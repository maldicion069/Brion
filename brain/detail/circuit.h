/* Copyright (c) 2013-2016, EPFL/Blue Brain Project
 *                          Juan Hernando <jhernando@fi.upm.es>
 *                          Adrien.Devresse@epfl.ch
 *                          Daniel.Nachbaur@epfl.ch
 *                          Stefan.Eilemann@epfl.ch
 *
 * This file is part of Brion <https://github.com/BlueBrain/Brion>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <brain/neuron/morphology.h>
#include <brion/blueConfig.h>
#include <brion/circuit.h>
#include <brion/morphology.h>
#include <brion/synapse.h>
#include <brion/synapseSummary.h>
#include <brion/target.h>
#include <brion/detail/lockHDF5.h>
#include <brion/detail/silenceHDF5.h>

#include <lunchbox/lock.h>
#include <lunchbox/lockable.h>
#include <lunchbox/log.h>
#include <lunchbox/persistentMap.h>
#include <lunchbox/scopedMutex.h>

#ifdef BRAIN_USE_MVD3
#  include <mvd/mvd3.hpp>
#  include <mvd/mvd_generic.hpp>
#endif

#include <boost/filesystem.hpp>
#include <boost/unordered_map.hpp>

#include <future>
#include <random>

namespace fs = boost::filesystem;
using boost::lexical_cast;

namespace brain
{

const std::string summaryFilename( "/nrn_summary.h5" );
const std::string afferentFilename( "/nrn.h5" );
const std::string efferentFilename( "/nrn_efferent.h5" );
const std::string afferentPositionsFilename( "/nrn_positions.h5" );
const std::string efferentPositionsFilename( "/nrn_positions_efferent.h5" );
const std::string extraFilename( "/nrn_extra.h5" );

namespace
{
#ifdef BRAIN_USE_MVD3
bool isSequence( const GIDSet& gids )
{
    return ( *gids.rbegin() - *gids.begin() + 1 ) == gids.size();
}

::MVD3::Range getRange( const GIDSet& gids )
{
    const size_t offset = ( *gids.begin( ));
    const size_t count = *gids.rbegin() - offset + 1;
    return ::MVD3::Range( offset - 1, count );
}

template< typename SrcArray, typename DstArray, typename AssignOp >
void assign( const ::MVD3::Range& range, const GIDSet& gids,
             SrcArray& src, DstArray& dst, const AssignOp& assignOp )
{
    if( isSequence( gids )) // OPT: no holes, no translation needed
    {
        std::transform( src.begin(), src.end(), dst.begin(), assignOp );
        return;
    }

    typename DstArray::iterator dst_it = dst.begin();
    for( GIDSet::const_iterator i = gids.begin(); i != gids.end(); ++i )
    {
        typename SrcArray::const_iterator src_it = src.begin();
        std::advance( src_it, *i - range.offset - 1 );
        *dst_it = assignOp( *src_it );
        ++dst_it;
    }
}

Vector3f toVector3f(
    const ::MVD3::Positions::const_subarray< 1 >::type& subarray )
{
    return Vector3f( subarray[0], subarray[1], subarray[2] );
}

Quaternionf toQuaternion(
    const ::MVD3::Rotations::const_subarray< 1 >::type& subarray )
{
    return Quaternionf( subarray[0], subarray[1], subarray[2], subarray[3] );
}

size_t nop( const size_t& in ) { return in; }
#endif

std::string toString( const std::string& in ) { return in; }
size_t toSize_t( const std::string& in ) { return std::stoul( in ); }

template<typename T>
void _shuffle( T& container )
{
    std::random_device randomDevice;
    std::mt19937_64 randomEngine( randomDevice( ));
    const char* seedEnv = getenv( "BRAIN_CIRCUIT_SEED" );
    if( seedEnv )
    {
        try
        {
            randomEngine.seed( std::stoul( seedEnv ));
        }
        catch( const std::exception& exc )
        {
            LBWARN << "Could not set BRAIN_CIRCUIT_SEED to " << seedEnv << ": "
                   << exc.what() << std::endl;
        }
    }
    std::shuffle( container.begin(), container.end(), randomEngine );
}

typedef boost::unordered_map< std::string, neuron::MorphologyPtr > Loaded;
} // anonymous namespace

class Circuit::Impl
{
public:
    explicit Impl( const brion::BlueConfig& config )
        : _circuitSource( config.getCircuitSource( ))
        , _morphologySource( config.getMorphologySource( ))
        , _synapseSource( config.getSynapseSource( ))
        , _targetSources( config.getTargetSources( ))
        , _cache( lunchbox::PersistentMap::createCache( ))
    {
    }

    virtual ~Impl() {}

    virtual size_t getNumNeurons() const = 0;

    const brion::URI& getCircuitSource() const
    {
        return _circuitSource;
    }

    GIDSet getGIDs() const
    {
        brain::GIDSet gids;
        brain::GIDSet::const_iterator hint = gids.begin();
        for( uint32_t i = 0; i < getNumNeurons(); ++i )
            hint = gids.insert( hint, i + 1 );
        return gids;
    }

    GIDSet getGIDs( const std::string& target ) const
    {
        if( _targetParsers.empty( ))
        {
            for( const URI& uri : _targetSources )
            {
                try
                {
                    _targetParsers.push_back( brion::Target( uri.getPath( )));
                }
                catch( const std::runtime_error& exc )
                {
                    LBWARN << "Failed to load targets from " << uri.getPath()
                           << ": " << exc.what() << std::endl;
                }
            }
        }
        return brion::Target::parse( _targetParsers, target );
    }

    GIDSet getRandomGIDs( const float fraction,
                          const std::string& target ) const
    {
        if( fraction < 0.f || fraction > 1.f )
            LBTHROW( std::runtime_error( "Fraction for getRandomGIDs() must be "
                                         "in the range [0,1]" ));

        const GIDSet& gids = target.empty() ? getGIDs() : getGIDs( target );
        uint32_ts randomGids( gids.begin(), gids.end( ));
        _shuffle( randomGids );
        randomGids.resize( size_t( std::ceil( randomGids.size() * fraction )));
        return GIDSet( randomGids.begin(), randomGids.end( ));
    }

    virtual Vector3fs getPositions( const GIDSet& gids ) const = 0;
    virtual size_ts getMTypes( const GIDSet& gids ) const = 0;
    virtual Strings getMorphologyNames() const = 0;
    virtual size_ts getETypes( const GIDSet& gids ) const = 0;
    virtual Strings getElectrophysiologyNames() const = 0;
    virtual Quaternionfs getRotations( const GIDSet& gids ) const = 0;
    virtual Strings getMorphologyNames( const GIDSet& gids ) const = 0;

    URI getMorphologyURI( const std::string& name ) const
    {
        URI uri;
        uri.setPath( _morphologySource.getPath() + "/" + name + ".h5" );
        uri.setScheme( "file" );
        return uri;
    }

    const brion::SynapseSummary& getSynapseSummary() const
    {
        lunchbox::ScopedWrite mutex( _synapseSummary );

        if( !(*_synapseSummary ))
            _synapseSummary->reset( new brion::SynapseSummary(
                               _synapseSource.getPath() + summaryFilename ));
        return **_synapseSummary;
    }

    const brion::Synapse& getSynapseAttributes( const bool afferent ) const
    {
        const size_t i = afferent ? 0 : 1;
        lunchbox::ScopedWrite mutex( _synapseAttributes[i] );

        if( !(*_synapseAttributes[i] ))
            _synapseAttributes[i]->reset(
                new brion::Synapse( _synapseSource.getPath() + (afferent ?
                                          afferentFilename : efferentFilename )));
        return **_synapseAttributes[i];
    }

    const brion::Synapse* getSynapseExtra() const
    {
        lunchbox::ScopedWrite mutex( _synapseExtra );

        if( !(*_synapseExtra ))
        {
            try
            {
                _synapseExtra->reset( new brion::Synapse(
                                 _synapseSource.getPath() + extraFilename ));
            }
            catch( ... )
            {
                return nullptr;
            }
        }
        return _synapseExtra->get();
    }

    const brion::Synapse& getSynapsePositions( const bool afferent ) const
    {
        const size_t i = afferent ? 0 : 1;
        lunchbox::ScopedWrite mutex( _synapsePositions[i] );

        if( !(*_synapsePositions[i] ))
            _synapsePositions[i]->reset(
                new brion::Synapse( _synapseSource.getPath() + (afferent ?
                      afferentPositionsFilename : efferentPositionsFilename )));
        return **_synapsePositions[i];
    }

    void saveToCache( const std::string& hash,
                      neuron::MorphologyPtr morphology ) const
    {
        if( _cache )
        {
            servus::Serializable::Data data = morphology->toBinary();
            _cache->insert( hash, data.ptr.get(), data.size );
        }
    }

    Loaded loadFromCache( const std::set< std::string >& hashes LB_UNUSED )
        const
    {
        Loaded loaded;
        if( _cache )
        {
            LBDEBUG << "Using cache for morphology loading" << std::endl;
            typedef std::future< std::pair< std::string,
                                            neuron::MorphologyPtr > > Future;
            std::vector< Future > futures;

            Strings keys( hashes.begin(), hashes.end( ));
            futures.reserve( keys.size( ));

            _cache->takeValues( keys, [&futures] ( const std::string& key,
                                                 char* data, const size_t size )
            {
                futures.push_back( std::async( std::launch::async,
                                               [key, data, size]
                {
                    neuron::MorphologyPtr morphology(
                                new neuron::Morphology( data, size ));
                    std::free( data );
                    return std::make_pair( key, morphology );
                }));
            });

            for( auto& future : futures )
                loaded.insert( future.get( ));

            LBINFO << "Loaded " << loaded.size() << " out of " << hashes.size()
                   << " morphologies from cache" << std::endl;
        }
        return loaded;
    }

    const brion::URI _circuitSource;
    const brion::URI _morphologySource;
    const brion::URI _synapseSource;
    const brion::URIs _targetSources;
    mutable brion::Targets _targetParsers;
    mutable lunchbox::PersistentMapPtr _cache;

    template< typename T >
    using LockPtr = lunchbox::Lockable< std::unique_ptr< T >>;

    mutable LockPtr< brion::SynapseSummary > _synapseSummary;
    mutable LockPtr< brion::Synapse > _synapseAttributes[2];
    mutable LockPtr< brion::Synapse > _synapseExtra;
    mutable LockPtr< brion::Synapse > _synapsePositions[2];
};

class MVD2 : public Circuit::Impl
{
public:
    MVD2( const brion::BlueConfig& config )
        : Impl( config )
        , _circuit( config.getCircuitSource().getPath( ))
    {}

    size_t getNumNeurons() const final
    {
        return _circuit.getNumNeurons();
    }

    Vector3fs getPositions( const GIDSet& gids ) const final
    {
        const brion::NeuronMatrix& data = _circuit.get( gids,
            brion::NEURON_POSITION_X | brion::NEURON_POSITION_Y |
            brion::NEURON_POSITION_Z );

        Vector3fs positions( gids.size( ));
#pragma omp parallel for
        for( size_t i = 0; i < gids.size(); ++i )
        {
            try
            {
                positions[i] =
                    brion::Vector3f( lexical_cast< float >( data[i][0] ),
                                     lexical_cast< float >( data[i][1] ),
                                     lexical_cast< float >( data[i][2] ));
            }
            catch( const boost::bad_lexical_cast& )
            {
                GIDSet::const_iterator gid = gids.begin();
                std::advance( gid, i );
                LBWARN << "Error parsing circuit position for gid "
                       << *gid << std::endl;
            }
        }
        return positions;
    }

    size_ts getMTypes( const GIDSet& gids ) const final
    {
        const brion::NeuronMatrix& matrix =  _circuit.get( gids,
                                                           brion::NEURON_MTYPE );
        size_ts result( matrix.shape()[ 0 ]);

        brion::NeuronMatrix::const_array_view<1>::type view =
            matrix[ boost::indices[brion::NeuronMatrix::index_range( )][ 0 ]];
        std::transform( view.begin(), view.end(), result.begin(), toSize_t );
        return result;
    }

    Strings getMorphologyNames() const final
    {
        return _circuit.getTypes( brion::NEURONCLASS_MTYPE );
    }

    size_ts getETypes( const GIDSet& gids ) const final
    {
        const brion::NeuronMatrix& matrix =  _circuit.get( gids,
                                                           brion::NEURON_ETYPE );
        size_ts result( matrix.shape()[ 0 ]);

        brion::NeuronMatrix::const_array_view<1>::type view =
            matrix[ boost::indices[brion::NeuronMatrix::index_range( )][ 0 ]];
        std::transform( view.begin(), view.end(), result.begin(), toSize_t );
        return result;
    }

    Strings getElectrophysiologyNames() const final
    {
        return _circuit.getTypes( brion::NEURONCLASS_ETYPE );
    }

    Quaternionfs getRotations( const GIDSet& gids ) const final
    {
        const float deg2rad = float( M_PI ) / 180.f;
        const brion::NeuronMatrix& data =
            _circuit.get( gids, brion::NEURON_ROTATION );
        Quaternionfs rotations( gids.size( ));

#pragma omp parallel for
        for( size_t i = 0; i < gids.size(); ++i )
        {
            try
            {
                // transform rotation Y angle in degree into rotation quaternion
                const float angle = lexical_cast<float>( data[i][0] ) * deg2rad;
                rotations[i] = Quaternionf( angle, Vector3f( 0, 1, 0 ));
            }
            catch( const boost::bad_lexical_cast& )
            {
                GIDSet::const_iterator gid = gids.begin();
                std::advance( gid, i );
                LBWARN << "Error parsing circuit orientation for gid "
                       << *gid << std::endl;
            }
        }
        return rotations;
    }

    Strings getMorphologyNames( const GIDSet& gids ) const final
    {
        const brion::NeuronMatrix& matrix =
                _circuit.get( gids, brion::NEURON_MORPHOLOGY_NAME );
        Strings result( matrix.shape()[ 0 ]);

        brion::NeuronMatrix::const_array_view<1>::type view =
            matrix[ boost::indices[brion::NeuronMatrix::index_range( )][ 0 ]];
        std::transform( view.begin(), view.end(), result.begin(), toString );
        return result;
    }

private:
    brion::Circuit _circuit;
};

#ifdef BRAIN_USE_MVD3
struct MVD3 : public Circuit::Impl
{
    MVD3( const brion::BlueConfig& config )
        : Impl( config )
        , _circuit( config.getCircuitSource().getPath( ))
    {}

    size_t getNumNeurons() const final
    {
        return _circuit.getNbNeuron();
    }

    Vector3fs getPositions( const GIDSet& gids ) const final
    {
        Vector3fs results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        try
        {
            brion::detail::SilenceHDF5 silence;
            lunchbox::ScopedWrite mutex( brion::detail::_hdf5Lock );
            const ::MVD3::Positions& positions = _circuit.getPositions( range );
            assign( range, gids, positions, results, toVector3f );
            return results;
        }
        catch( const HighFive::Exception& e )
        {
            LBTHROW( std::runtime_error( "Exception in getPositions(): " +
                                         std::string( e.what( ))));
        }
    }

    size_ts getMTypes( const GIDSet& gids ) const final
    {
        size_ts results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        try
        {
            brion::detail::SilenceHDF5 silence;
            lunchbox::ScopedWrite mutex( brion::detail::_hdf5Lock );
            const size_ts& mtypes = _circuit.getIndexMtypes( range );
            assign( range, gids, mtypes, results, nop );
            return results;
        }
        catch( const HighFive::Exception& e )
        {
            LBTHROW( std::runtime_error( "Exception in getMTypes(): " +
                                         std::string( e.what( ))));
        }
    }

    Strings getMorphologyNames() const final
    {
        return _circuit.listAllMtypes();
    }

    size_ts getETypes( const GIDSet& gids ) const final
    {
        size_ts results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        try
        {
            brion::detail::SilenceHDF5 silence;
            lunchbox::ScopedWrite mutex( brion::detail::_hdf5Lock );
            const size_ts& etypes = _circuit.getIndexEtypes( range );
            assign( range, gids, etypes, results, nop );
            return results;
        }
        catch( const HighFive::Exception& e )
        {
            LBTHROW( std::runtime_error( "Exception in getETypes(): " +
                                         std::string( e.what( ))));
        }
    }

    Strings getElectrophysiologyNames() const final
    {
        return _circuit.listAllEtypes();
    }

    Quaternionfs getRotations( const GIDSet& gids ) const final
    {
        Quaternionfs results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        try
        {
            brion::detail::SilenceHDF5 silence;
            lunchbox::ScopedWrite mutex( brion::detail::_hdf5Lock );
            const ::MVD3::Rotations& rotations = _circuit.getRotations( range );
            assign( range, gids, rotations, results, toQuaternion );
            return results;
        }
        catch( const HighFive::Exception& e )
        {
            LBTHROW( std::runtime_error( "Exception in getRotations(): " +
                                         std::string( e.what( ))));
        }
    }

    Strings getMorphologyNames( const GIDSet& gids ) const final
    {
        Strings results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        try
        {
            brion::detail::SilenceHDF5 silence;
            lunchbox::ScopedWrite mutex( brion::detail::_hdf5Lock );
            const Strings& morphos = _circuit.getMorphologies( range );
            assign( range, gids, morphos, results, toString );
            return results;
        }
        catch( const HighFive::Exception& e )
        {
            LBTHROW( std::runtime_error( "Exception in getMorphologyNames(): " +
                                         std::string( e.what( ))));
        }
    }

private:
    ::MVD3::MVD3File _circuit;
};
#endif

}
