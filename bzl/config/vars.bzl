load("@bazel_skylib//lib:selects.bzl", "selects")

# ate_pairing
# CFLAGS_WARN = ["-Wall", "-Wextra", "-Wformat=2",
#                "-Wcast-qual", "-Wcast-align", "-Wwrite-strings",
#                "-Wfloat-equal", "-Wpointer-arith"]

WARNINGS = ["-Wall", "-Wextra", "-Wno-unused-variable"]

DEBUG_FLAGS = select({
    "//:enable_debug": ["-g"],
    "//conditions:default": ["-g0"]
}) + select({
    "//:with_verbose": ["-v"],
    "//conditions:default": []
}) + select({
    "//:macos_disable_debug": ["-UDEBUG"],
    "//conditions:default": []
})

OPTIMIZE_CXXFLAGS = select({
    "@//bzl/config:enable_optimization": ["-flto", "-fuse-linker-plugin", "-O2"],
    "//conditions:default": ["-O0"]
})

OPTIMIZE_LINKFLAGS = select({
    "@//bzl/config:enable_optimization": ["-flto", "-fuse-linker-plugin"],
    "//conditions:default": []
})

CPPFLAGS = DEBUG_FLAGS + WARNINGS
CFLAGS   = []
CXXFLAGS = OPTIMIZE_CXXFLAGS
LDFLAGS  = []
#######################
####    DEFINES    ####
DDEBUG = select({
    "//:enable_debug": ["DEBUG"],
    "//conditions:default": [] # ["NDEBUG"]
})

DBINARY_OUTPUT = select({
    "@//bzl/config:enable_binary_output": ["BINARY_OUTPUT"],
    "//conditions:default": []
})

DCURVE = select({
    "@//:curve_bn128": ["CURVE_BN128"],
    "@//:curve_alt_bn128": ["CURVE_ALT_BN128"],
    "@//:curve_edwards": ["CURVE_EDWARDS"],
    "@//:curve_mnt4": ["CURVE_MNT4"],
    "@//:curve_mnt6": ["CURVE_MNT6"],
    "//conditions:default": ["CURVE_BN128"]
})

# transitive, for ate_pairing:
DLIBGMP = select({
    "@//bzl/config:enable_libgmp" : ["MIE_ATE_USE_GMP"],
    "//conditions:default": []
})

DLOWMEM = select({
    "@//bzl/config:enable_lowmem": ["LOWMEM"],
    "//conditions:default": []
})

DMONTGOMERY_OUTPUT = select({
    "@//bzl/config:enable_montgomery_output": ["MONTGOMERY_OUTPUT"],
    "//conditions:default": []
})

DMULTICORE = select({
    "@//bzl/config:enable_multicore": ["MULTICORE"],
    "//conditions:default": []
})

DNO_PT_COMPRESSION = select({
    "@//bzl/config:disable_point_compression": ["NO_PT_COMPRESSION"],
    "//conditions:default": []
})

DPROFILE_OP_COUNTS = select({
    "@//bzl/config:enable_profile_op_counts": ["PROFILE_OP_COUNTS"],
    "//conditions:default": []
})

DUSE_MIXED_ADDITION = select({
    "@//bzl/config:enable_mixed_addition": ["USE_MIXED_ADDITION"],
    "//conditions:default": []
})

DNO_PROCPS = selects.with_or({
    ("@//bzl/config:disable_procps", "@//bzl/host:macos"): ["NO_PROCPS"],
    "//conditions:default": []
})

DCXX_DEBUG = select({
    "@//bzl/config:enable_cxx_debug": ["_GLIBCXX_DEBUG", "_GLIBCXX_DEBUG_PEDANTIC"],
    "//conditions:default": []
})

DOPTIMIZE = select({
    "@//bzl/config:enable_optimization": ["NDEBUG"],
    "//conditions:default": []
})

DUSE_ASM = select({
    "@//bzl/config:enable_asm": ["USE_ASM"],
    "//conditions:default": []
})
