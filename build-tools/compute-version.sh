#!/bin/bash
#
# Create PEP-440 compliant version identifiers using information stored in
# git.
#
# Uses the following logic:
#
# * If we're not in a git repository (or git is not installed), return the
#   contents of the file named `VERSION`.  If VERSION does not exist, return
#   "unknown".
#
# * If the current commit has an annotated tag attached to it, the version
#   is simply the tag with the leading 'v' removed.
#
# * If the current commit is past an annotated tag, the version is constructed at:
#
#    '{tag}.post{commitcount}+{gitsha}'
#
#  where {commitcount} is the number of commits after the tag (obtained with `git describe`)
#
# * If there are no annotated tags in the past, the version is:
#
#    '1.0.2.post{commitcount}+{gitsha}'
#
#  where 1.0.2 is picked up from DEFAULT_VERSION variable.
#
# * If the working directory is dirty, add '.dirty' to the version.
#
# * If '-u' is given on the command line, update the file VERSION if the
#   computed version differs from the one currently in the file.
#
# Inspired by https://github.com/pyfidelity/setuptools-git-version
#
# Author: mjuric@astro.washington.edu (http://github.com/mjuric)
#

DEFAULT_VERSION=v1.0.2

compute_version() {
	# See if we have any tags defined
	GD=$(git describe --always --long --dirty --match 'v[0-9]*' 2>/dev/null)
	if [[ ! -z $GD ]]; then
		IFS=- read TAG COUNT HASH DIRTY <<< "$GD"

		# No tags at all
		if [[ -z "$HASH" ]]; then
			# Number of commits since the beginning
			HASH="g$TAG"
			DIRTY="$COUNT"
			TAG="$DEFAULT_VERSION"
			COUNT="$(git rev-list --count HEAD)"
		fi

		# remove the 'v' prefix from tag
		TAG="${TAG#v}"

		# remove the 'g' prefix from hash
		HASH="${HASH#g}"

		# handle the dirty flag
		if [[ ! -z $DIRTY ]]; then
			# If dirty, return $HASH.dirty, or just .dirty (if $HASH is empty)
			[[ ! -z "HASH" ]] && HASH="$HASH."
			HASH="$HASH$DIRTY"
		fi

		# construct the version
		VERSION="$TAG"
		if [[ $COUNT != "0" ]]; then
			VERSION="$VERSION.post$COUNT"
		fi
	
		if [[ $COUNT != "0" || $DIRTY ]]; then
			VERSION="$VERSION+$HASH"
		fi

		echo "$VERSION"
	elif [[ -f VERSION ]]; then
		# If a file named "VERSION" exists, use its contents as the version
		head -n 1 VERSION
	else
		echo "unknown"
	fi
}

# Assuming this script is in $ROOT/build-tools/ subdirectory;
# Change the directory to $ROOT.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR/.."

# Print out the version
VERSION="$(compute_version)"
echo "$VERSION"

# Save the version to a file, if requested (and if different from what's already in the file)
if [[ $1 == "-u" ]]; then
	if [[ ! -f VERSION || $VERSION != $(cat VERSION) ]]; then
		echo "$VERSION" > VERSION
	fi
fi
