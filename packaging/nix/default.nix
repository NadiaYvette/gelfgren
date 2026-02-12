{ lib
, stdenv
, rustPlatform
, cbindgen
, pkg-config
}:

rustPlatform.buildRustPackage rec {
  pname = "gelfgren";
  version = "0.1.0";

  src = ../..;

  cargoLock = {
    lockFile = ../../Cargo.lock;
  };

  nativeBuildInputs = [
    cbindgen
    pkg-config
  ];

  postBuild = ''
    cd gelfgren-ffi
    cbindgen --config cbindgen.toml --output $out/include/gelfgren.h
    cd ..
  '';

  postInstall = ''
    mkdir -p $out/lib/pkgconfig
    cat > $out/lib/pkgconfig/gelfgren.pc << EOF
prefix=$out
exec_prefix=\''${prefix}
libdir=\''${exec_prefix}/lib
includedir=\''${prefix}/include

Name: Gelfgren
Description: Piecewise rational interpolation library
Version: ${version}
Libs: -L\''${libdir} -lgelfgren
Cflags: -I\''${includedir}
EOF
  '';

  meta = with lib; {
    description = "Piecewise rational interpolation library";
    longDescription = ''
      Gelfgren is a high-performance numerical computing library implementing
      piecewise rational interpolation methods based on Jan Gelfgren's 1975
      research. Written in Rust with bindings for 17 programming languages.
    '';
    homepage = "https://github.com/yourusername/gelfgren";
    license = with licenses; [ mit asl20 ];
    maintainers = [ "Nadia Chambers <nadia.chambers@iohk.io>" ];
    platforms = platforms.unix;
  };
}
