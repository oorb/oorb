{
  description = "oorb: an open source orbit computation package";

  inputs.nixpkgs.url     = "github:nixos/nixpkgs/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.flake-compat    = { url = "github:edolstra/flake-compat"; flake = false; };

  outputs = { self, nixpkgs, flake-compat, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkg-name = "oorb";
        pkgs = nixpkgs.legacyPackages.${system};
          
        eg = pkgs.runCommand
          "oorb"
          { preferLocalBuild = true; buildInputs = [ pkg-name ]; }
          '''';

        revision = "${self.lastModifiedDate}-${self.shortRev or "dirty"}";
      in {
        defaultPackage = self.packages.${system}.${pkg-name};

        devShell = pkgs.mkShell {
          buildInputs = [
            pkgs.curl
            pkgs.fswatch
            pkgs.gcc13
            pkgs.gfortran13
            pkgs.gnuplot
            pkgs.lapack
            pkgs.python3
            pkgs.texlive.combined.scheme-basic
            pkgs.pkg-config
            pkgs.sourceHighlight
            pkgs.stdenv.cc.cc.lib
            pkgs.zlib
            (pkgs.python3.withPackages (ps: [
              ps.numpy
              ps.pyp
              ps.pytest
            ]))
          ];

          shellHook = ''
            export SHELL=/run/current-system/sw/bin/bash
            export LANG=en_US.UTF-8
            export PS1="nix|$PS1"
            export OORB_DATA=./data
          '';
        };
      }
    );
}
