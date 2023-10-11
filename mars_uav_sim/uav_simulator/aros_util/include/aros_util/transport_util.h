#ifndef TRANSPORT_IMAGE_H
#define TRANSPORT_IMAGE_H

#include "opencv2/opencv.hpp"
#include "sensor_msgs/Image.h"
#include "sensor_msgs/fill_image.h"
#include "sensor_msgs/image_encodings.h"
#include <boost/regex.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <ros/package.h>
#include <ros/ros.h>

/**
 * copy depthStrToInt() function from cv_bridge.cpp
 */
static int
depthStrToInt( const std::string depth )
{
    if ( depth == "8U" )
    {
        return 0;
    }
    else if ( depth == "8S" )
    {
        return 1;
    }
    else if ( depth == "16U" )
    {
        return 2;
    }
    else if ( depth == "16S" )
    {
        return 3;
    }
    else if ( depth == "32S" )
    {
        return 4;
    }
    else if ( depth == "32F" )
    {
        return 5;
    }
    return 6;
}

/**
 * copy getCvType() function from cv_bridge.cpp
 * @param encoding
 * @return
 */
static int
getCvType( const std::string& encoding )
{
    // Check for the most common encodings first
    if ( encoding == sensor_msgs::image_encodings::BGR8 )
        return CV_8UC3;
    if ( encoding == sensor_msgs::image_encodings::MONO8 )
        return CV_8UC1;
    if ( encoding == sensor_msgs::image_encodings::RGB8 )
        return CV_8UC3;
    if ( encoding == sensor_msgs::image_encodings::MONO16 )
        return CV_16UC1;
    if ( encoding == sensor_msgs::image_encodings::BGR16 )
        return CV_16UC3;
    if ( encoding == sensor_msgs::image_encodings::RGB16 )
        return CV_16UC3;
    if ( encoding == sensor_msgs::image_encodings::BGRA8 )
        return CV_8UC4;
    if ( encoding == sensor_msgs::image_encodings::RGBA8 )
        return CV_8UC4;
    if ( encoding == sensor_msgs::image_encodings::BGRA16 )
        return CV_16UC4;
    if ( encoding == sensor_msgs::image_encodings::RGBA16 )
        return CV_16UC4;

    // For bayer, return one-channel
    if ( encoding == sensor_msgs::image_encodings::BAYER_RGGB8 )
        return CV_8UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_BGGR8 )
        return CV_8UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_GBRG8 )
        return CV_8UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_GRBG8 )
        return CV_8UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_RGGB16 )
        return CV_16UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_BGGR16 )
        return CV_16UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_GBRG16 )
        return CV_16UC1;
    if ( encoding == sensor_msgs::image_encodings::BAYER_GRBG16 )
        return CV_16UC1;

    // Miscellaneous
    if ( encoding == sensor_msgs::image_encodings::YUV422 )
        return CV_8UC2;

    // Check all the generic content encodings
    boost::cmatch m;

    if ( boost::regex_match( encoding.c_str( ), m, boost::regex( "(8U|8S|16U|16S|32S|32F|64F)C([0-9]+)" ) ) )
    {
        return CV_MAKETYPE( depthStrToInt( m[1].str( ) ), atoi( m[2].str( ).c_str( ) ) );
    }

    if ( boost::regex_match( encoding.c_str( ), m, boost::regex( "(8U|8S|16U|16S|32S|32F|64F)" ) ) )
    {
        return CV_MAKETYPE( depthStrToInt( m[1].str( ) ), 1 );
    }

    throw std::invalid_argument( "Unrecognized image encoding [" + encoding + "]" );
}
/**
 * Reference: http://stackoverflow.com/questions/4239993/determining-endianness-at-compile-time
 * TODO: This might be not safe!!!
 * @return
 */
static bool
is_little_endian( )
{
    short int n = 0x1;
    return ( *( char* )&n == 1 );
}
static bool
is_big_endian( )
{
    union {
        uint32_t i;
        char c[4];
    } bint = { 0x01020304 };
    // cout << bool(bint.c[0] == 1) << "==============" << endl;
    return ( bint.c[0] == 1 );
}
/**
 * TODO: copy matFromImage() function from cv_bridge.cpp, not fully the same due to boost::endian
 * @param source
 * @return
 */

static void
toImageMsg( sensor_msgs::Image& ros_image, cv::Mat image )
{
    //  ros_image.header = header;
    ros_image.height = image.rows;
    ros_image.width  = image.cols;
    //  ros_image.encoding = encoding;
    //  ros_image.is_bigendian = (boost::endian::order::native == boost::endian::order::big);
    ros_image.is_bigendian = is_big_endian( );
    ros_image.step         = image.cols * image.elemSize( );
    size_t size            = ros_image.step * image.rows;
    ros_image.data.resize( size );

    if ( image.isContinuous( ) )
    {
        memcpy( ( char* )( &ros_image.data[0] ), image.data, size );
    }
    else
    {
        // Copy by row by row
        uchar* ros_data_ptr = ( uchar* )( &ros_image.data[0] );
        uchar* cv_data_ptr  = image.data;
        for ( int i = 0; i < image.rows; ++i )
        {
            memcpy( ros_data_ptr, cv_data_ptr, ros_image.step );
            ros_data_ptr += ros_image.step;
            cv_data_ptr += image.step;
        }
    }
}

static cv::Mat
cvmat_from_msg( const std::string& encoding, const uint32_t step, const uint8_t is_bigendian, const uint32_t height, const uint32_t width, const uchar* data )
{
    int source_type  = getCvType( encoding );
    int byte_depth   = sensor_msgs::image_encodings::bitDepth( encoding ) / 8;
    int num_channels = sensor_msgs::image_encodings::numChannels( encoding );

    // If the endianness is the same as locally, share the data
    cv::Mat mat( height, width, source_type, const_cast<uchar*>( data ), step );

    if ( ( !is_little_endian( ) && is_bigendian ) || ( is_little_endian( ) && !is_bigendian ) || ( byte_depth == 1 ) )
        return mat;

    // Otherwise, reinterpret the data as bytes and switch the channels accordingly
    mat = cv::Mat( height, width, CV_MAKETYPE( CV_8U, num_channels * byte_depth ), const_cast<uchar*>( data ), step );
    cv::Mat mat_swap( height, width, mat.type( ) );

    std::vector<int> fromTo;
    fromTo.reserve( num_channels * byte_depth );
    for ( int i = 0; i < num_channels; ++i )
        for ( int j = 0; j < byte_depth; ++j )
        {
            fromTo.push_back( byte_depth * i + j );
            fromTo.push_back( byte_depth * i + byte_depth - 1 - j );
        }
    cv::mixChannels( std::vector<cv::Mat>( 1, mat ), std::vector<cv::Mat>( 1, mat_swap ), fromTo );

    // Interpret mat_swap back as the proper type
    mat_swap = cv::Mat( height, width, source_type, mat_swap.data, mat_swap.step );
    return mat_swap;
}

static cv::Mat
matFromImage( const sensor_msgs::Image& source )
{
    return cvmat_from_msg( source.encoding, source.step, source.is_bigendian, source.height, source.width, &source.data[0] );
}

static cv::Mat
matFromImage( const sensor_msgs::ImageConstPtr& source )
{
    return cvmat_from_msg( source->encoding, source->step, source->is_bigendian, source->height, source->width, &source->data[0] );
}

#endif // TRANSPORT_IMAGE_H
