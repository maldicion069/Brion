
/* Copyright (c) 2013-2016, EPFL/Blue Brain Project
 *                          bbp-open-source@googlegroups.com
 *                          Daniel Nachbaur <daniel.nachbaur@epfl.ch>
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

#include "synapse.h"
#include "detail/lockHDF5.h"
#include "detail/silenceHDF5.h"

#include <bitset>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/unordered_map.hpp>

#include <H5Cpp.h>

#include <lunchbox/atomic.h>
#include <lunchbox/log.h>
#include <lunchbox/persistentMap.h>
#include <lunchbox/scopedMutex.h>

#include <fstream>


namespace brion
{
namespace detail
{
namespace
{
static lunchbox::a_ssize_t _cacheHits;
static lunchbox::a_ssize_t _cacheMiss;
}

struct Dataset
{
    H5::DataSet dataset;
    H5::DataSpace dataspace;
    hsize_t dims[2];
};
namespace fs = boost::filesystem;
using boost::lexical_cast;

/** Access a single synapse file (nrn*.h5 or nrn*.h5.<int> */
class SynapseFile : public boost::noncopyable
{
public:
    explicit SynapseFile( const std::string& source )
        : _cache( lunchbox::PersistentMap::createCache( ))
        , _cacheKey( fs::canonical( fs::path( source )).generic_string( ))
    {
        lunchbox::ScopedWrite mutex( detail::_hdf5Lock );

        try
        {
            SilenceHDF5 silence;
            _file.openFile( source, H5F_ACC_RDONLY );
        }
        catch( const H5::Exception& exc )
        {
            LBTHROW( std::runtime_error( "Could not open synapse file " + source
                                         + ": " + exc.getDetailMsg( )));
        }

        Dataset dataset;
        const std::string& datasetName = _file.getObjnameByIdx( 0 );
        if( !_openDataset( datasetName, dataset ))
            LBTHROW( std::runtime_error( "Cannot open dataset in synapse file "
                                         + source ));

        _numAttributes = dataset.dims[1];
        if( _numAttributes != SYNAPSE_ALL &&
            _numAttributes != SYNAPSE_POSITION_ALL &&
            _numAttributes != 1 /* nrn_extra */)
        {
            LBTHROW( std::runtime_error( source + " not a valid synapse file"));
        }
    }

    ~SynapseFile()
    {
        lunchbox::ScopedWrite mutex( detail::_hdf5Lock );

        _file.close();
    }

    template< size_t N >
    SynapseMatrix read( const uint32_t gid, const uint32_t attributes ) const
    {
        const std::bitset< N > bits( attributes );
        if( !bits.any( ))
            return SynapseMatrix();

        std::string cacheKey;
        if( _cache )
        {
            lunchbox::ScopedWrite mutex( _cacheLock );
            cacheKey = _cacheKey + "/" + lexical_cast< std::string >( gid ) +
                       "/" + lexical_cast< std::string >( attributes );
            const std::string& cached = (*_cache)[ cacheKey ];
            if( !cached.empty( ))
            {
                if( (++_cacheHits % 5000) == 0 )
                    LBDEBUG << int( float( _cacheHits ) /
                                    float( _cacheHits+_cacheMiss )*100.f + .5f )
                            << "% cache hit rate" << std::endl;

                const size_t dim0 = cached.size() / bits.count() /
                                    sizeof( float );
                SynapseMatrix values( boost::extents[ dim0 ][ bits.count( )]);
                ::memcpy( values.data(), cached.data(), cached.size( ));
                return values;
            }
            if( (++_cacheMiss % 5000) == 0 )
                LBDEBUG << int( float( _cacheHits ) /
                                float( _cacheHits + _cacheMiss ) * 100.f + .5f )
                        << "% cache hit rate" << std::endl;
        }

        lunchbox::ScopedWrite mutex( detail::_hdf5Lock );
        Dataset dataset;
        if( !_openDataset( gid, dataset ))
            return SynapseMatrix();

        dataset.dataspace.selectNone();
        for( size_t i = 0; i < bits.size(); ++i )
        {
            if( bits.test( i ))
            {
                const hsize_t readCounts[2] = { dataset.dims[0], 1 };
                const hsize_t readOffsets[2] = { 0, i };
                dataset.dataspace.selectHyperslab( H5S_SELECT_OR, readCounts,
                                                   readOffsets );
            }
        }

        SynapseMatrix values( boost::extents[dataset.dims[0]][bits.count()] );
        const hsize_t targetSizes[2] = { dataset.dims[0], bits.count() };
        H5::DataSpace targetspace( 2, targetSizes );

        dataset.dataset.read( values.data(), H5::PredType::NATIVE_FLOAT,
                              targetspace, dataset.dataspace );

        if( _cache )
        {
            lunchbox::ScopedWrite cacheMutex( _cacheLock );
            const size_t size = dataset.dims[0] * bits.count() * sizeof( float );
            if( !_cache->insert( cacheKey, values.data(), size ))
                LBWARN << "Failed to insert synapse information for GID " << gid
                       << " into cache; item size is " << float(size) / LB_1MB
                       << " MB" << std::endl;

        }
        return values;
    }

    size_t getNumSynapses( const GIDSet& gids ) const
    {
        lunchbox::ScopedWrite mutex( detail::_hdf5Lock );

        size_t numSynapses = 0;
        BOOST_FOREACH( const uint32_t gid, gids )
        {
            Dataset dataset;
            if( !_openDataset( gid, dataset ))
                continue;
            numSynapses += dataset.dims[0];
        }
        return numSynapses;
    }

    bool _openDataset( const uint32_t gid, Dataset& dataset ) const
    {
        std::stringstream name;
        name << "a" << gid;
        return _openDataset( name.str(), dataset );
    }

    bool _openDataset( const std::string& name, Dataset& dataset ) const
    {
        try
        {
            SilenceHDF5 silence;
            dataset.dataset = _file.openDataSet( name );
        }
        catch( const H5::Exception& )
        {
            LBVERB << "Could not find synapse dataset for " << name << ": "
                   << std::endl;
            return false;
        }

        dataset.dataspace = dataset.dataset.getSpace();
        if( dataset.dataspace.getSimpleExtentNdims() != 2 )
        {
            LBERROR << "Synapse dataset is not 2 dimensional" << std::endl;
            return false;
        }

        if( dataset.dataspace.getSimpleExtentDims( dataset.dims ) < 0 )
        {
            LBERROR << "Synapse dataset dimensions could not be retrieved"
                    << std::endl;
            return false;
        }

        return true;
    }

    SynapseMatrix read( const uint32_t gid, const uint32_t attributes ) const
    {
        switch( _numAttributes )
        {
        case SYNAPSE_ALL:
            return read< SYNAPSE_ALL >( gid, attributes );
        case SYNAPSE_POSITION_ALL:
            return read< SYNAPSE_POSITION_ALL >( gid, attributes );
        case 1:
            // nrn_extra
            return read< 1 >( gid, 1 );
        default:
            LBERROR << "Synapse file " << _file.getFileName()
                    << " has unknown number of attributes: " << _numAttributes
                    << std::endl;
            return SynapseMatrix();
        }
    }

private:
    lunchbox::PersistentMapPtr _cache;
    mutable lunchbox::Lock _cacheLock;
    std::string _cacheKey;
    H5::H5File _file;
    size_t _numAttributes;
};

/** Implement the logic to read a merged .h5 or individual .h5.<int> files */
class Synapse : public boost::noncopyable
{
public:
    explicit Synapse( const std::string& source )
        : _file( 0 )
        , _gid( 0 )
    {
        try
        {
            _file = new SynapseFile( source );
        }
        catch( const std::runtime_error& )
        {
            LBINFO << "No merged synapse file found at " << source << std::endl;

            const fs::path sourcePath( source );
            const fs::path dir = sourcePath.parent_path();
            const std::string filename = sourcePath.filename().generic_string();

            _createIndex( dir, filename );
            _fillFilemap( dir, filename );
        }
    }

    ~Synapse()
    {
        delete _file;
    }

    SynapseMatrix read( const uint32_t gid, const uint32_t attributes ) const
    {
        if( _findFile( gid ))
            return _file->read( gid, attributes );
        return SynapseMatrix();
    }

    size_t getNumSynapses( const GIDSet& gids ) const
    {
        size_t numSynapses = 0;
        BOOST_FOREACH( const uint32_t gid, gids )
        {
            if( !_findFile( gid ))
                continue;

            GIDSet set;
            set.insert( gid );
            numSynapses += _file->getNumSynapses( set );
        }
        return numSynapses;
    }

private:
    typedef boost::unordered_map< uint32_t, std::string > GidFileMap;
    typedef boost::unordered_map< std::string, H5::H5File > UnmappedFiles;

    mutable SynapseFile* _file;
    mutable uint32_t _gid; // current or 0 for all
    mutable UnmappedFiles _unmappedFiles;
    mutable GidFileMap _fileMap;

    void _createIndex( const fs::path& dir, const std::string& filename )
    {
        // extract the GID->file mapping from the merge_nrn.sh script
        const bool afferent = filename.find( "efferent" ) == std::string::npos;
        const fs::path merge_nrn = dir / (afferent ? "merge_nrn.sh"
                                                   : "merge_nrn_efferent.sh");
        const std::ifstream mergeFile( merge_nrn.generic_string().c_str( ));
        if( !mergeFile.is_open( ))
        {
            LBWARN << "No merge file found in " << dir
                   << " to build lookup index; loading data will be slow"
                   << std::endl;
            return;
        }

        std::stringstream buffer;
        buffer << mergeFile.rdbuf();

        const boost::regex commentregx( "#.*?\\n" );
        const std::string content = boost::regex_replace( buffer.str(),
                                                          commentregx , "\n" );

        const boost::regex regx( "\\$CMD -i \\$H5.(?<number>[0-9]+) -o "
                                 "\\$H5 -s /a(?<gid>[0-9]+)" );
        const int subs[] = {1, 2};
        boost::sregex_token_iterator i( content.begin(), content.end(), regx,
                                        subs );

        const std::string basename =
                            (dir / fs::path( filename )).generic_string() + ".";
        for( boost::sregex_token_iterator end; i != end; )
        {
            const std::string& fileNumber = *i++;
            const uint32_t fileGID = boost::lexical_cast< uint32_t >(*i++);
            _fileMap[ fileGID ] = basename + fileNumber;
        }
    }

    void _fillFilemap( const fs::path& dir, const std::string& filename )
    {
        // try to open in individual files
        const boost::regex filter( filename + "\\.[0-9]+$" );
        fs::directory_iterator end;
        for( fs::directory_iterator i( dir ); i != end; ++i )
        {
            const fs::path candidate = i->path().filename();
            boost::smatch match;

            if( !boost::filesystem::is_regular_file( i->status( )) ||
                !boost::regex_match( candidate.string(), match, filter ))
            {
                continue;
            }

            _unmappedFiles.insert( std::make_pair( i->path().string(),
                                                   H5::H5File( )));
        }

        if( _unmappedFiles.empty( ))
        {
            LBTHROW( std::runtime_error( "Could not find synapse files " +
                                         dir.string() + "/" + filename ));
        }
    }

    bool _findFile( const uint32_t gid ) const
    {
        if( _file && ( _gid == gid || _gid == 0 ))
            return true;

        const std::string& filename = _findFilename( gid );
        if( filename.empty( ))
            return false;

        delete _file;
        _file = new SynapseFile( filename );
        _gid = gid;
        return true;
    }

    std::string _findFilename( const uint32_t gid ) const
    {
        if( !_fileMap[ gid ].empty( ))
            return _fileMap[ gid ];

        // at this point we can only search in each file for the GID which
        // usually results in waiting for I/O and non-parallizable search thanks
        // to HDF5

        lunchbox::ScopedWrite mutex( detail::_hdf5Lock );
        SilenceHDF5 silence;

        BOOST_FOREACH( UnmappedFiles::value_type& entry, _unmappedFiles )
        {
            const std::string& candidate = entry.first;
            H5::H5File& file = entry.second;
            // keeping the files open 'saves' some time
            if( file.getId() <= 0 )
                file.openFile( candidate, H5F_ACC_RDONLY );

            try
            {
                std::stringstream name;
                name << "a" << gid;

                // this trial-and-error is the 'fastest' path found
                file.openDataSet( name.str( ));
                _fileMap[ gid ] = candidate;
                return candidate;
            }
            catch( const H5::Exception& ) {}
        }
        return std::string();
    }
};

}

Synapse::Synapse( const std::string& source )
    : _impl( new detail::Synapse( source ))
{
}

Synapse::~Synapse()
{
    delete _impl;
}

SynapseMatrix Synapse::read( const uint32_t gid,
                             const uint32_t attributes ) const
{
    return _impl->read( gid, attributes );
}

size_t Synapse::getNumSynapses( const GIDSet& gids ) const
{
    return _impl->getNumSynapses( gids );
}

}
