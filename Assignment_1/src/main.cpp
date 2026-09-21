////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    double det_vector = u[0] * v[1] - u[1] * v[0];
    return det_vector;
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    double det1 = det(b - a, c - a);
    double det2 = det(b - a, d - a);
    double det3 = det(d - c, a - c);
    double det4 = det(d - c, b - c);

    auto on_segment = [](const Vector2d &p, const Vector2d &u, const Vector2d &v)
    {
        return p.x() >= std::min(u.x(), v.x()) && p.x() <= std::max(u.x(), v.x()) &&
               p.y() >= std::min(u.y(), v.y()) && p.y() <= std::max(u.y(), v.y());
    };

    if (det1 == 0 && on_segment(c, a, b)) return true;
    if (det2 == 0 && on_segment(d, a, b)) return true;
    if (det3 == 0 && on_segment(a, c, d)) return true;
    if (det4 == 0 && on_segment(b, c, d)) return true;

    return ((det1 < 0 && det2 > 0) || (det1 > 0 && det2 < 0)) &&
           ((det3 < 0 && det4 > 0) || (det3 > 0 && det4 < 0));
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    if (poly.empty())
    {
        return false;
    }
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    double min_x = poly[0][0];
    double max_x = poly[0][0];
    double min_y = poly[0][1];
    double max_y = poly[0][1];

    for (size_t i = 1; i < poly.size(); ++i)
    {
        min_x = std::min(min_x, poly[i][0]);
        max_x = std::max(max_x, poly[i][0]);
        min_y = std::min(min_y, poly[i][1]);
        max_y = std::max(max_y, poly[i][1]);
    }

    Vector2d outside(max_x + (max_x - min_x) + 1.0, max_y + (max_y - min_y) + 1.0);

    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    if (query.x() < min_x || query.x() > max_x || query.y() < min_y || query.y() > max_y)
    {
        return false;
    }

    int intersection_count = 0;

    for (size_t i = 0; i < poly.size(); ++i)
    {
        const Vector2d &a = poly[i];
        const Vector2d &b = poly[(i + 1) % poly.size()];

        if (intersect_segment(a, b, query, query))
            return true;

        if ((a.y() > query.y()) != (b.y() > query.y()))
        {
            double x = a.x() + (query.y() - a.y()) *
                                (b.x() - a.x()) / (b.y() - a.y());

            if (x > query.x() && x <= outside.x())
                ++intersection_count;
        }
    
    }
    return (intersection_count % 2 == 1);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);

    int count;

    if (!(in >> count) || count < 0)
    {
        throw std::runtime_error("Invalid XYZ header from: " + filename);
    }

    for (int i = 0; i < count; ++i)
    {
        double x, y, z;
        if (!(in >> x >> y >> z))
            throw std::runtime_error("Invalid XYZ point from: " + filename);

        points.push_back(Vector2d(x, y));
    }
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open output: " << filename << '\n';
        return;
    }
    out << points.size() << '\n';
    for (size_t i = 0; i < points.size(); ++i)
    {
        out << points[i][0] << " " << points[i][1] << " " << 0.0 << std::endl;
    }
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}
