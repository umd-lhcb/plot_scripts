let
  pkgs = import <nixpkgs> {};
in

pkgs.mkShell {
  name = "plot_scripts";
  buildInputs = with pkgs; [
    stdenv
    python2
    root
  ];
}
