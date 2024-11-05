#!/bin/bash
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
cd "$SCRIPT_DIR" || { echo "Failed to change directory to $SCRIPT_DIR"; exit 1; }

check_install(){
	if perl -M"$1" -e '1;' 2>/dev/null; then
		echo "Package $1 is installed."
	else
		cpan $1
fi
}

# install perl pkgs
check_install "List::MoreUtils"
check_install "Getopt::Long"
check_install "MIME::Base64"
check_install "Statistics::TTest"
check_install "Text::NSP::Measures::2D::Fisher::twotailed"
check_install "Statistics::Multtest"
check_install "File::Which"

# install deAPA
R -e "install.packages('./deAPA_1.0.tar.gz',repos = NULL, type = 'source')"