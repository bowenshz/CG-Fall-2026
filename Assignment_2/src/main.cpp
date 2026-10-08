// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);
    const double sphere_radius = 0.9;
    const Vector3d sphere_center(0, 0, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;
            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            // TODO change this with the generic case
            // Vector2d ray_on_xy(ray_origin(0), ray_origin(1));

            const Vector3d origin_offset = ray_origin - sphere_center;
            const double direction_squared = ray_direction.squaredNorm();
            const double projected_offset = ray_direction.dot(origin_offset);
            const double discriminant = projected_offset * projected_offset - direction_squared * (origin_offset.squaredNorm() - sphere_radius * sphere_radius);
            
            if (direction_squared == 0 || discriminant < 0)
            {
                continue;
            }
            
            double intersection_distance = (-projected_offset - sqrt(discriminant)) / direction_squared;
            
            if (intersection_distance < 0)
            {
                intersection_distance = (-projected_offset + sqrt(discriminant)) / direction_squared;
            }

            // if (ray_on_xy.norm() < sphere_radius)
            
            if (intersection_distance >= 0)
            {
                // The ray hit the sphere, compute the exact intersection point
                // Vector3d ray_intersection(
                //     ray_on_xy(0), ray_on_xy(1),
                //     sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));
                
                Vector3d ray_intersection = ray_origin + intersection_distance * ray_direction;

                // Compute normal at the intersection point
                // Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            
            Matrix3d intersection_matrix;
            intersection_matrix.col(0) = pgram_u;
            intersection_matrix.col(1) = pgram_v;
            intersection_matrix.col(2) = -ray_direction;
            
            const FullPivLU<Matrix3d> intersection_solver(intersection_matrix);
            
            if (!intersection_solver.isInvertible())
            {
                continue;
            }
            const Vector3d intersection_parameters = intersection_solver.solve(ray_origin - pgram_origin);
            // if (true)
            if (intersection_parameters(0) >= 0 && intersection_parameters(0) <= 1 && intersection_parameters(1) >= 0 && intersection_parameters(1) <= 1 && intersection_parameters(2) >= 0)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection
                // point
                // Vector3d ray_intersection(0, 0, 0);
                Vector3d ray_intersection = ray_origin + intersection_parameters(2) * ray_direction;

                // TODO: Compute normal at the intersection point
                // Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }
    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(1, 0.4, 0);
    const Vector3d pgram_v(0, 0.7, -10);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            // const Vector3d ray_origin = pixel_center;
            // const Vector3d ray_direction = camera_view_direction;
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = (pixel_center - camera_origin).normalized();

            // TODO: Check if the ray intersects with the parallelogram
            
            Matrix3d intersection_matrix;
            intersection_matrix.col(0) = pgram_u;
            intersection_matrix.col(1) = pgram_v;
            intersection_matrix.col(2) = -ray_direction;
            
            const FullPivLU<Matrix3d> intersection_solver(intersection_matrix);
            
            if (!intersection_solver.isInvertible())
            {
                continue;
            }
            
            const Vector3d intersection_parameters = intersection_solver.solve(ray_origin - pgram_origin);
            
            // if (true)
            if (intersection_parameters(0) >= 0 && intersection_parameters(0) <= 1 && intersection_parameters(1) >= 0 && intersection_parameters(1) <= 1 && intersection_parameters(2) >= 0)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection point
                // Vector3d ray_intersection(0, 0, 0);
                Vector3d ray_intersection = ray_origin + intersection_parameters(2) * ray_direction;

                // TODO: Compute normal at the intersection point
                // Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = pgram_u.cross(pgram_v).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }
    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    // MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd R = MatrixXd::Zero(800, 800);
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    const Vector3d light_intesity(1, 1, 1);
    double ambient = 0.1;

    // for (unsigned i = 0; i < C.cols(); ++i)
    for (unsigned i = 0; i < A.cols(); ++i)
    {
        // for (unsigned j = 0; j < C.rows(); ++j)
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            // const Vector3d ray_origin = pixel_center;
            // const Vector3d ray_direction = camera_view_direction;
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = (pixel_center - camera_origin).normalized();

            // Intersect with the sphere
            // TODO: implement the generic ray sphere intersection
            const Vector3d origin_offset = ray_origin - sphere_center;
            const double direction_squared = ray_direction.squaredNorm();
            const double projected_offset = ray_direction.dot(origin_offset);
            const double discriminant = projected_offset * projected_offset - direction_squared * (origin_offset.squaredNorm() - sphere_radius * sphere_radius);
            
            if (direction_squared == 0 || discriminant < 0)
            {
                continue;
            }
            
            double intersection_distance = (-projected_offset - sqrt(discriminant)) / direction_squared;
            
            if (intersection_distance < 0)
            {
                intersection_distance = (-projected_offset + sqrt(discriminant)) / direction_squared;
            }
            
            // if (true)
            
            if (intersection_distance >= 0)
            {
                // TODO: The ray hit the sphere, compute the exact intersection point
                // Vector3d ray_intersection(0, 0, 0);
                Vector3d ray_intersection = ray_origin + intersection_distance * ray_direction;

                // TODO: Compute normal at the intersection point
                // Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                // TODO: Add shading parameter here
                // const double diffuse = (light_position - ray_intersection).normalized().dot(ray_normal);
                // const double specular = (light_position - ray_intersection).normalized().dot(ray_normal);
                const Vector3d light_direction = (light_position - ray_intersection).normalized();
                const Vector3d view_direction = -ray_direction;
                const Vector3d halfway_direction = (light_direction + view_direction).normalized();
                const double diffuse = std::max(light_direction.dot(ray_normal), 0.);
                const double specular = diffuse > 0 ? pow(std::max(halfway_direction.dot(ray_normal), 0.), specular_exponent) : 0.;

                // Simple diffuse model
                // C(i, j) = ambient + diffuse + specular;
                const Vector3d pixel_color = ambient * diffuse_color +
                    (diffuse * diffuse_color + specular * specular_color).cwiseProduct(light_intesity);

                // Clamp to zero
                // C(i, j) = std::max(C(i, j), 0.);
                R(i, j) = std::max(pixel_color(0), 0.);
                G(i, j) = std::max(pixel_color(1), 0.);
                B(i, j) = std::max(pixel_color(2), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }
    // Save to png
    // write_matrix_to_png(C, C, C, A, filename);
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
