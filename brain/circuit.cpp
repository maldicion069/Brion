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

#ifdef _WINDOWS
// first thing to include to get M_PI defined
#define _USE_MATH_DEFINES
#include <cmath>
#endif

#include "circuit.h"
#include "detail/circuit.h"

#include "synapsesStream.h"

namespace brain
{

<<<<<<< HEAD
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
#endif

std::string toString( const std::string& in ) { return in; }
}

class Circuit::Impl
{
public:
    explicit Impl( const brion::BlueConfig& config )
        : _morphologySource( config.getMorphologySource( ))
        , _targetSources( config.getTargetSources( ))
    {}
    virtual ~Impl() {}

    virtual size_t getNumNeurons() const = 0;

    GIDSet getGIDs() const
    {
        brain::GIDSet gids;
        brain::GIDSet::const_iterator hint = gids.begin();
        for( size_t i = 0; i < getNumNeurons(); ++i )
            hint = gids.insert( hint, uint32_t( i ) + 1 );
        return gids;
    }

    GIDSet getGIDs( const std::string& target ) const
    {
        if( _targetParsers.empty( ))
        {
            BOOST_FOREACH( const URI& uri, _targetSources )
                _targetParsers.push_back( brion::Target( uri.getPath( )));
        }
        return brion::Target::parse( _targetParsers, target );
    }

    virtual Vector3fs getPositions( const GIDSet& gids ) const = 0;
    virtual Quaternionfs getRotations( const GIDSet& gids ) const = 0;
    virtual Strings getMorphologyNames( const GIDSet& gids ) const = 0;

    URI getMorphologyURI( const std::string& name ) const
    {
        URI uri;
        uri.setPath( _morphologySource.getPath() + "/" + name + ".h5" );
        uri.setScheme( "file" );
        return uri;
    }

private:
    const brion::URI _morphologySource;
    const brion::URIs _targetSources;
    mutable brion::Targets _targetParsers;
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
        for( int i = 0; i < int( gids.size()); ++i )
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

    Quaternionfs getRotations( const GIDSet& gids ) const final
    {
        const float deg2rad = float( M_PI ) / 180.f;
        const brion::NeuronMatrix& data =
            _circuit.get( gids, brion::NEURON_ROTATION );
        Quaternionfs rotations( gids.size( ));

#pragma omp parallel for
        for( int i = 0; i < int( gids.size()); ++i )
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
        brion::NeuronMatrix matrix =
            _circuit.get( gids, brion::NEURON_MORPHOLOGY_NAME );
        Strings res( matrix.shape()[ 0 ]);

        brion::NeuronMatrix::array_view<1>::type view =
            matrix[ boost::indices[brion::NeuronMatrix::index_range( )][ 0 ]];
        std::transform( view.begin(), view.end(), res.begin(), toString );
        return res;
    }

private:
    brion::Circuit _circuit;
};

#ifdef BRAIN_USE_MVD3
class MVD3 : public Circuit::Impl
{
public:
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
        const ::MVD3::Positions& positions = _circuit.getPositions( range );
        assign( range, gids, positions, results, toVector3f );
        return results;
    }

    Quaternionfs getRotations( const GIDSet& gids ) const final
    {
        Quaternionfs results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        const ::MVD3::Rotations& rotations = _circuit.getRotations( range );
        assign( range, gids, rotations, results, toQuaternion );
        return results;
    }

    Strings getMorphologyNames( const GIDSet& gids ) const final
    {
        Strings results( gids.size( ));
        const ::MVD3::Range& range = getRange( gids );
        const Strings& morphos = _circuit.getMorphologies( range );
        assign( range, gids, morphos, results, toString );
        return results;
    }

private:
    ::MVD3::MVD3File _circuit;
};
#endif

=======
>>>>>>> eb47d43e3626a519b128628c3def5e9088a4b184
Circuit::Impl* newImpl( const brion::BlueConfig& config )
{
    const std::string circuit = config.getCircuitSource().getPath();
    if( boost::algorithm::ends_with( circuit, ".mvd2" ))
        return new MVD2( config );
#ifdef BRAIN_USE_MVD3
    return new MVD3( config );
#else
    throw std::runtime_error( "MVD3 support requires CMake 3" );
#endif
}

Circuit::Circuit( const URI& source )
  : _impl( newImpl( brion::BlueConfig( source.getPath( ))))
{
}

Circuit::Circuit( const brion::BlueConfig& config )
  : _impl( newImpl( config ))
{
}

Circuit::~Circuit()
{
}

GIDSet Circuit::getGIDs() const
{
    return _impl->getGIDs();
}

GIDSet Circuit::getGIDs( const std::string& target ) const
{
    return _impl->getGIDs( target );
}

GIDSet Circuit::getRandomGIDs( const float fraction ) const
{
    return _impl->getRandomGIDs( fraction, "" );
}

GIDSet Circuit::getRandomGIDs( const float fraction,
                               const std::string& target ) const
{
    return _impl->getRandomGIDs( fraction, target );
}

URIs Circuit::getMorphologyURIs( const GIDSet& gids ) const
{
    const Strings& names = _impl->getMorphologyNames( gids );

    URIs uris;
    uris.reserve( names.size( ));
    for( Strings::const_iterator i = names.begin(); i < names.end(); ++i )
       uris.push_back( _impl->getMorphologyURI( *i ));
    return uris;
}

neuron::Morphologies Circuit::loadMorphologies( const GIDSet& gids,
                                                const Coordinates coords ) const
{
    const URIs& uris = getMorphologyURIs( gids );

    // < GID, hash >
    Strings gidHashes;
    gidHashes.reserve( uris.size( ));
    std::set< std::string > hashes;
    GIDSet::const_iterator gid = gids.begin();
    for( size_t i = 0; i < uris.size(); ++i, ++gid )
    {
        std::string hash( fs::canonical( uris[i].getPath( )).generic_string( ));

        if( coords == Coordinates::global )
        {
            // store circuit + GID for transformed morphology
            hash += fs::canonical(
                        _impl->getCircuitSource().getPath( )).generic_string() +
                        boost::lexical_cast< std::string >( *gid );
            hash = servus::make_uint128( hash ).getString();
        }
        else
            hash = servus::make_uint128( hash ).getString();

        gidHashes.push_back( hash );
        hashes.insert( hash );
    }

    Loaded loaded = _impl->loadFromCache( hashes );

    // resolve missing morphologies and put them in GID-order into result
    neuron::Morphologies result;
    result.reserve( uris.size( ));

    const Matrix4fs transforms =
             coords == Coordinates::global ? getTransforms( gids ) : Matrix4fs();
    for( size_t i = 0; i < uris.size(); ++i )
    {
        const URI& uri = uris[i];

        const std::string& hash = gidHashes[i];
        Loaded::const_iterator it = loaded.find( hash );
        if( it == loaded.end( ))
        {
            neuron::MorphologyPtr morphology;
            const brion::Morphology raw( uri.getPath( ));
            if( coords == Coordinates::global )
                morphology.reset( new neuron::Morphology( raw, transforms[i] ));
            else
                morphology.reset( new neuron::Morphology( raw ));

            loaded.insert( std::make_pair( hash, morphology ));

            _impl->saveToCache( hash, morphology );

            result.push_back( morphology );
        }
        else
            result.push_back( it->second );
    }

    return result;
}

Vector3fs Circuit::getPositions( const GIDSet& gids ) const
{
    return _impl->getPositions( gids );
}

size_ts Circuit::getMorphologyTypes( const GIDSet& gids ) const
{
    return _impl->getMTypes( gids );
}

Strings Circuit::getMorphologyNames() const
{
    return _impl->getMorphologyNames();
}

size_ts Circuit::getElectrophysiologyTypes( const GIDSet& gids ) const
{
    return _impl->getETypes( gids );
}

Strings Circuit::getElectrophysiologyNames() const
{
    return _impl->getElectrophysiologyNames();
}

Matrix4fs Circuit::getTransforms( const GIDSet& gids ) const
{
    const Vector3fs& positions = _impl->getPositions( gids );
    const Quaternionfs& rotations = _impl->getRotations( gids );
    if( positions.size() != rotations.size( ))
        throw std::runtime_error(
            "Positions not equal rotations for given GIDs" );

    Matrix4fs transforms( positions.size( ));

#pragma omp parallel for
    for( int i = 0; i < int( positions.size()); ++i )
        transforms[i] = Matrix4f( rotations[i], positions[i] );
    return transforms;

}

Quaternionfs Circuit::getRotations( const GIDSet& gids ) const
{
    return _impl->getRotations( gids );
}

size_t Circuit::getNumNeurons() const
{
    return _impl->getNumNeurons();
}

SynapsesStream Circuit::getAfferentSynapses( const GIDSet& gids,
                                          const SynapsePrefetch prefetch ) const
{
    return SynapsesStream( *this, gids, true, prefetch );
}

SynapsesStream Circuit::getEfferentSynapses( const GIDSet& gids,
                                          const SynapsePrefetch prefetch ) const
{
    return SynapsesStream( *this, gids, false, prefetch );
}

SynapsesStream Circuit::getProjectedSynapses( const GIDSet& preGIDs,
                                              const GIDSet& postGIDs,
                                          const SynapsePrefetch prefetch ) const
{
    return SynapsesStream( *this, preGIDs, postGIDs, prefetch );
}

}
