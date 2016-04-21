# EupsPkg config file. Sourced by 'eupspkg'

_ensure_exists()
{
	hash "$1" 2>/dev/null || die "need '$1' to install this product. please install it and try again."
}

prep()
{
	# check for system prerequisites
	_ensure_exists wget
	_ensure_exists gfortran

	default_prep
}

config()
{
	./configure gfortran opt
}

build()
{
	PKGROOT="$PWD"

	# built & test OpenOrb
	( cd main && make oorb )

	# update JPL Ephemeris files and make 405 and 430 ephemeris files
        ( export EPH_TYPE=405 )
	( cd data/JPL_ephemeris && make && make test )
        ( export EPH_TYPE=430 )
	( cd data/JPL_ephemeris && make && make test )

	# build & test python bindings
	(
		cd python
		make pyoorb

		export OORB_CONF="$PKGROOT/main"
		export OORB_DATA="$PKGROOT/data"
		export LD_LIBRARY_PATH="$PKGROOT/lib:$LD_LIBRARY_PATH"
		export DYLD_LIBRARY_PATH="$PKGROOT/lib:$DYLD_LIBRARY_PATH"
		python test.py
	)

}

install()
{
	clean_old_install
	
	make PREFIX="$PREFIX" install

	install_ups
}
