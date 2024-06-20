using GLMakie

"""
    Cone(quality::Number)

Return mesh::Mesh{3, Float32, Triangle}
Implementation of a cone mesh/marker using GLMakie primitive functions using quality setting (sampling).
Adapted (and accordingly edited) from `Beautiful Makie` example file: https://beautiful.makie.org/examples/dashboards/matcap
"""
Cone(q::Number) = merge([
    Makie._circle(Point3f(0), 0.5f0, Vec3f(0, 0, -1), q),
    Makie._mantle(Point3f(0), Point3f(0, 0, 1), 0.5f0, 0.0f0, q)])
