function vec_angle(v1, v2)
    n1 = norm(v1)
    n2 = norm(v2)
    return 2atan(norm(n1*v2 - n2*v1), norm(n1*v2 + n2*v1))
end

const I3 = SMatrix{3, 3}((
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    ))

xrot(x) = SMatrix{3, 3}((
    1.0, 0.0, 0.0,
    0.0, cos(x), sin(x),
    0.0, -sin(x), cos(x)
    ))

zrot(x) = SMatrix{3, 3}((
    cos(x), sin(x), 0.0,
    -sin(x), cos(x), 0.0,
    0.0, 0.0, 1.0,
    ))