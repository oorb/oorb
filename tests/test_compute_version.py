##
## Unit tests for compute-version.sh script (to be run with pytest)
##

from __future__ import print_function
from subprocess import check_output
import os, os.path, tempfile, shutil
from pkg_resources import parse_version

def shell(cmd):
	return check_output(cmd, shell=True).decode('utf-8')

def check_version(verexp):
	#
	# Hand-parse the versions
	#
	verall = check_output("./build-tools/compute-version.sh", shell=True).decode('utf-8').rstrip()

	# SHAs at the end will be changing (they depend on dates, author, etc..). Split them off
	(verexp, shaexp) = verexp.split('+') if '+' in verexp else (verexp, None)
	(ver, sha)       = verall.split('+') if '+' in verall else (verall, None)

	# Split sha on '.dirty', if exists
	(shaexp, dirtyexp) = shaexp.split('.') if shaexp is not None and '.' in shaexp else (shaexp, None)
	(sha, dirty)       = sha.split('.')    if sha    is not None and '.' in sha    else (sha, None)

	#
	# Check that different parts match expectations
	#
	# Version strings have to match, and the (non)existance of SHAs has to be the same
	assert ver == verexp and (sha is not None) == (shaexp is not None)

	# The dirty flags have to match, and spell 'dirty'
	assert dirty == dirtyexp
	assert dirty == None or dirty == 'dirty'

	# If we have a sha, it has to be the same as the current commit SHA
	if sha is not None:
		shacur = shell("git rev-parse HEAD")[:len(sha)]
		print((sha, shacur))
		assert sha == shacur

	return parse_version(verall)

def test_compute_version():
	cwd = os.getcwd()
	try:
		# Create a mock git repository in a temporary dir,
		# copy the tested script there.
		gitwd = tempfile.mkdtemp()
		btdir = os.path.join(gitwd, 'build-tools')
		os.mkdir(btdir)
		shutil.copy2('build-tools/compute-version.sh', btdir)
		os.chdir(gitwd)

		#
		# Running outside of a git repository, with no VERSION file
		#
		ver = check_version("unknown")

		#
		# Running outside of a git repository, with a VERSION file
		#
		shell("""
			echo 1.0.0 > VERSION
		""")
		ver, prev = check_version("1.0.0"), ver
		assert ver > prev

		#
		# initialize a git repo, with a single commit
		#
		shell("""
			git init
			echo A dummy test repo >> README.md
			git add README.md
			git commit -m 'Initial commit'
		""")
		ver, prev = check_version("1.0.2.post1+a1310fc"), ver
		assert ver > prev

		#
		# check the 'dirty' flag
		#
		shell("""
			echo "FOO" >> README.md
		""")
		ver, prev = check_version("1.0.2.post1+a1310fc.dirty"), ver

		#
		# Now with two commits
		#
		shell("""
			echo "A dummy file A" >> A
			git add A
			git commit -m 'A' -a
		""")
		ver, prev = check_version("1.0.2.post2+53dca7c"), ver
		assert ver > prev

		#
		# make sure we're NOT picking up non-annotated tags
		#
		shell("git tag v1.5.0")
		ver, prev = check_version("1.0.2.post2+53dca7c"), ver
		assert ver == prev

		#
		# make sure we're NOT picking up tags that don't begin
		# with a 'v'
		shell("git tag 1.6.0 -a -m 'Version 1.6.0'")
		ver, prev = check_version("1.0.2.post2+53dca7c"), ver
		assert ver == prev

		#
		# Run on annotated tag, properly formatted
		#
		shell("git tag v2.0.0 -a -m 'Version v2.0.0'")
		ver, prev = check_version("2.0.0"), ver
		assert ver > prev

		#
		# now with one additional commit
		#
		shell("echo foo >> A; git commit -m 'added to A' -a")
		ver, prev = check_version("2.0.0.post1+8e51bcc"), ver
		assert ver > prev

		#
		# add one more commit, to ensure the version increments
		#
		shell("echo foo >> A; git commit -m 'added to A' -a")
		ver, prev = check_version("2.0.0.post2+0b24346"), ver
		assert ver > prev
	finally:
		os.chdir(cwd)
