#!/bin/bash

# Determine if the current commit is tagged
tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0

# Install Latex to build the documentation
sudo apt-get install --no-install-recommends \
  texlive-base texlive-latex-base texlive-generic-recommended \
  texlive-fonts-recommended texlive-fonts-extra texlive-extra-utils \
  texlive-latex-recommended texlive-latex-extra texinfo lmodern

# If it is, build the .tar.gz
echo "Building release $tag..."
make release || exit 1
release=`ls jaatha_*.tar.gz` || exit 1

# Configure git and make the repo use https-auth
echo "Setting up git..."
git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"
git remote rm origin
git remote add origin "https://${GH_TOKEN}@github.com/paulstaab/jaatha"
git fetch origin releases
git checkout releases

# Add the release 
git add $release

# And push everything to GitHub
echo "Pushing to GitHub..."
git commit -m "Add release $tag"
git push -q origin releases 2>&1 > /dev/null || exit 1
echo "Done"
