!host_build|!cross_compile {
    QMAKE_CXXFLAGS+=-F/opt/aarch64-apple-darwin20/aarch64-apple-darwin20/sys-root/System/Library/Frameworks
    QMAKE_RANLIB=/opt/bin/aarch64-apple-darwin20-libgfortran5-cxx11/ranlib
    QMAKE_MACOSX_DEPLOYMENT_TARGET=11.0
}
host_build {
    QT_CPU_FEATURES.x86_64 = mmx sse sse2
} else {
    QT_CPU_FEATURES.arm64 = neon
}
QT.global_private.enabled_features = alloca_h alloca dbus dlopen gui network reduce_exports relocatable sql system-zlib testlib widgets xml
QT.global_private.disabled_features = sse2 alloca_malloc_h android-style-assets avx2 private_tests dbus-linked gc_binaries intelcet libudev posix_fallocate reduce_relocations release_tools stack-protector-strong zstd
QMAKE_LIBS_LIBDL = 
QT_COORD_TYPE = double
QMAKE_LIBS_ZLIB = -lz
CONFIG += cross_compile compile_examples largefile neon precompile_header
QT_BUILD_PARTS += libs
QT_HOST_CFLAGS_DBUS += 
