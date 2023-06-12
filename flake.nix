{
  description = "eg haskell code";

  inputs.nixpkgs.url     = "github:nixos/nixpkgs/nixpkgs-23.05";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.flake-compat    = { url = "github:edolstra/flake-compat"; flake = false; };

  outputs = { self, nixpkgs, flake-compat, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkg-name = "oorb";
        pkgs = nixpkgs.legacyPackages.${system};
        haskell-pkgs = pkgs.haskell.packages.ghc961;

        eg = pkgs.runCommand
          "eg"
          { preferLocalBuild = true; buildInputs = [ pkg-name ]; }
          '''';

        revision = "${self.lastModifiedDate}-${self.shortRev or "dirty"}";
        jailbreak-unbreak = pkg: pkgs.haskell.lib.doJailbreak (pkg.overrideAttrs (_: { meta = { }; }));
      in {
        defaultPackage = self.packages.${system}.${pkg-name};

        devShell = pkgs.mkShell {
          buildInputs = [
            pkgs.pkg-config
            pkgs.sourceHighlight
            pkgs.stdenv.cc.cc.lib
            pkgs.zlib
          ];

          inputsFrom = builtins.attrValues self.packages.${system};

          shellHook = ''
            export LANG=en_US.UTF-8
            export PS1="nix|-$PS1"
          '';
        };
      }
    );
}
