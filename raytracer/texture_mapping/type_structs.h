//
// Created by furkan on 12/31/25.
//

#ifndef TYPE_STRUCTS_H
#define TYPE_STRUCTS_H

enum class TextureType : int {
    image = 0,
    perlin = 1,
    checkerboard = 2
};

enum class Interpolation : int {
    nearest,
    bilinear,
    trilinear
};

enum class DecalMode : int {
    replace_kd,
    blend_kd,
    replace_ks,
    replace_background,
    replace_normal,
    bump_normal,
    replace_all
};

enum class NoiseConversion : bool {
    linear,
    absval
};

#endif //TYPE_STRUCTS_H
