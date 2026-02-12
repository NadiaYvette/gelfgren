{
  description = "Gelfgren piecewise rational interpolation library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, rust-overlay }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        overlays = [ (import rust-overlay) ];
        pkgs = import nixpkgs {
          inherit system overlays;
        };

        gelfgren = pkgs.callPackage ./default.nix { };

      in
      {
        packages.default = gelfgren;
        packages.gelfgren = gelfgren;

        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            rust-bin.stable.latest.default
            cbindgen
            pkg-config

            # Language-specific tools
            ocaml
            opam
            julia
            go
            mlton
            sbcl
            guile
            swi-prolog
          ];

          shellHook = ''
            export LD_LIBRARY_PATH=${gelfgren}/lib:$LD_LIBRARY_PATH
            echo "Gelfgren development environment"
            echo "Library: ${gelfgren}/lib"
            echo "Headers: ${gelfgren}/include"
          '';
        };

        apps.default = flake-utils.lib.mkApp {
          drv = gelfgren;
        };
      }
    );
}
