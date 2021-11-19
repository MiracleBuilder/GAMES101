/**
 * @file rasterizer.cpp
 * @author your name (you@domain.com)
 * @brief comment the raterizer function of draw, is no msaa
 * @version 0.1
 * @date 2021-11-19
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>

rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f &v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

// added by user
static bool isSameSide(const Eigen::Vector3f& pt1, const Eigen::Vector3f& pt2, 
                       const Eigen::Vector3f& a, const Eigen::Vector3f& b)
{
    // reference: https://blackpawn.com/texts/pointinpoly/default.html
    Eigen::Vector3f cp1 = (b-a).cross(pt1-a);
    Eigen::Vector3f cp2 = (b-a).cross(pt2-a);
    if (cp1.dot(cp2)>0)
    {
        return true;
    }
    else
    {
        return false;
    }    
}

// function interface has been modified for msaa bonus
// static bool insideTriangle(int x, int y, const Vector3f* _v)
static bool insideTriangleSameSide(float x, float y, const Vector3f* _v)
{   
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    // reference: https://blackpawn.com/texts/pointinpoly/default.html
    Eigen::Vector3f pt;
    pt << (float)x, (float)y, 0.0f;
    bool isInside = isSameSide(pt, _v[0], _v[1], _v[2]) && isSameSide(pt, _v[1], _v[0], _v[2]) && isSameSide(pt, _v[2], _v[0], _v[1]);
    return isInside;
}

static bool insideTriangle(int x, int y, const Vector3f *_v)
{
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    // _v 0 1 2对应a b c
    Eigen::Vector3f oo(x + 0.5, y + 0.5, 0);
    Eigen::Vector3f ab = _v[1] - _v[0];
    Eigen::Vector3f ao = oo - _v[0];
    Eigen::Vector3f ac = _v[2] - _v[0];

    Eigen::Vector3f bc = _v[2] - _v[1];
    Eigen::Vector3f bo = oo - _v[1];

    Eigen::Vector3f ca = _v[0] - _v[2];
    Eigen::Vector3f co = oo - _v[2];

    // 同号即可
    if (ab.cross(ao).z() > 0 && bc.cross(bo).z() > 0 && ca.cross(co).z() > 0 ||
        ab.cross(ao).z() < 0 && bc.cross(bo).z() < 0 && ca.cross(co).z() < 0)
    {
        return true;
    }

    return false;
}

static bool insideTriangleFloat(float x, float y, const Vector3f *_v)
{
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    // _v 0 1 2对应a b c
    Eigen::Vector3f oo(x, y, 0);
    Eigen::Vector3f ab = _v[1] - _v[0];
    Eigen::Vector3f ao = oo - _v[0];
    Eigen::Vector3f ac = _v[2] - _v[0];

    Eigen::Vector3f bc = _v[2] - _v[1];
    Eigen::Vector3f bo = oo - _v[1];

    Eigen::Vector3f ca = _v[0] - _v[2];
    Eigen::Vector3f co = oo - _v[2];

    // 同号即可
    if (ab.cross(ao).z() > 0 && bc.cross(bo).z() > 0 && ca.cross(co).z() > 0 ||
        ab.cross(ao).z() < 0 && bc.cross(bo).z() < 0 && ca.cross(co).z() < 0)
    {
        return true;
    }

    return false;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f *v)
{
    float c1 = (x * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * y + v[1].x() * v[2].y() - v[2].x() * v[1].y()) / (v[0].x() * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * v[0].y() + v[1].x() * v[2].y() - v[2].x() * v[1].y());
    float c2 = (x * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * y + v[2].x() * v[0].y() - v[0].x() * v[2].y()) / (v[1].x() * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * v[1].y() + v[2].x() * v[0].y() - v[0].x() * v[2].y());
    float c3 = (x * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * y + v[0].x() * v[1].y() - v[1].x() * v[0].y()) / (v[2].x() * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * v[2].y() + v[0].x() * v[1].y() - v[1].x() * v[0].y());
    return {c1, c2, c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto &buf = pos_buf[pos_buffer.pos_id];
    auto &ind = ind_buf[ind_buffer.ind_id];
    auto &col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto &i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
            mvp * to_vec4(buf[i[0]], 1.0f),
            mvp * to_vec4(buf[i[1]], 1.0f),
            mvp * to_vec4(buf[i[2]], 1.0f)};
        //Homogeneous division
        for (auto &vec : v)
        {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto &vert : v)
        {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            // z的值要取反
            vert.z() = -vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        // rasterize_triangle(t);
        rasterize_triangle_msaa(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle &t)
{
    // TODO : Find out the bounding box of current triangle.
    auto v = t.toVector4();
    // the range of x
    int boxLeft = std::floor(std::min(v[0].x(), std::min(v[1].x(), v[2].x())));
    int boxRight = std::ceil(std::max(v[0].x(), std::max(v[1].x(), v[2].x())));

    // the range of y
    int boxTop = std::ceil(std::max(v[0].y(), std::max(v[1].y(), v[2].y())));
    int boxBottom = std::floor(std::min(v[0].y(), std::min(v[1].y(), v[2].y())));

    // iterate through the pixel and find if the current pixel is inside the triangle
    for (int x = boxLeft; x <= boxRight; x++)
    {
        for (int y = boxBottom; y <= boxTop; y++)
        {
            if (insideTriangle(x, y, t.v))
            {
                // get the interpolated z value
                float alpha, beta, gamma;
                // what's the meaning of computeBarycentric2D?
                std::tie(alpha, beta, gamma) = computeBarycentric2D(x, y, t.v);
                float wReciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float zInterpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                zInterpolated *= wReciprocal;

                // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted
                if (zInterpolated < depth_buf[get_index(x, y)])
                {
                    set_pixel(Vector3f(x, y, zInterpolated), t.getColor());
                    depth_buf[get_index(x, y)] = zInterpolated;
                }
            }
        }
    }
}

// the step of msaa sampling
float stepX[4] = {0.25, 0.75, 0.25, 0.75};
float stepY[4] = {0.25, 0.25, 0.75, 0.75};

void rst::rasterizer::rasterize_triangle_msaa(const Triangle &t)
{
    // TODO : Find out the bounding box of current triangle.
    auto v = t.toVector4();
    // the range of x
    int boxLeft = std::floor(std::min(v[0].x(), std::min(v[1].x(), v[2].x())));
    int boxRight = std::ceil(std::max(v[0].x(), std::max(v[1].x(), v[2].x())));

    // the range of y
    int boxTop = std::ceil(std::max(v[0].y(), std::max(v[1].y(), v[2].y())));
    int boxBottom = std::floor(std::min(v[0].y(), std::min(v[1].y(), v[2].y())));

    // iterate through the pixel and find if the current pixel is inside the triangle
    for (int x = boxLeft; x <= boxRight; x++)
    {
        for (int y = boxBottom; y <= boxTop; y++)
        {
            int count = 0;
            for (int i : {0, 1})
            {
                for (int j : {0, 1})
                {
                    float newX = x + stepX[i];
                    float newY = y + stepY[j];
                    if (insideTriangleFloat(newX, newY, t.v))
                    {
                        // get the interpolated z value
                        float alpha, beta, gamma;
                        // what's the meaning of computeBarycentric2D?
                        std::tie(alpha, beta, gamma) = computeBarycentric2D(newX, newY, t.v);
                        float wReciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                        float zInterpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                        zInterpolated *= wReciprocal;

                        int buffIndex = x * 2 + i + (y * 2 + j) * width * 2;
                        if (zInterpolated < depth_buf_msaa_2x2[buffIndex])
                        {
                            count++;
                            depth_buf_msaa_2x2[buffIndex] = zInterpolated;
                        }
                    }
                }
            }
            if (count > 0)
            {
                float percentage = count / 4.0f;
                set_pixel(Vector3f(x, y, 1.0f), t.getColor() * percentage);
            }
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f &m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f &v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f &p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
        // msaa
        std::fill(depth_buf_msaa_2x2.begin(), depth_buf_msaa_2x2.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
    // the depth of msaa initial
    depth_buf_msaa_2x2.resize(w * 2 * h * 2);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height - 1 - y) * width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f &point, const Eigen::Vector3f &color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;
}

void rst::rasterizer::add_pixel(const Eigen::Vector3f &point, const Eigen::Vector3f &color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] += color;
}

// clang-format on