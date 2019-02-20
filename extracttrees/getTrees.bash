#!/bin/bash
if [[ $# -ne 1 ]]
then
    echo "usage: bash $0 <ntree>"
else
    ntree=$1
    # generate trees from scikit learn predictor files
    # this write C++ code necessary to import the predictor into MoleTools-based code
    # MUST READ: https://en.wikipedia.org/wiki/Tree_traversal#Depth-first_search_2

    python printTrees.py > data/trees

    echo "//Aaron T. Frank"
    echo "//Sean M. Law"
    echo ""
    echo "/*"
    echo "This file is part of MoleTools."
    echo ""
    echo "MoleTools is free software: you can redistribute it and/or modify"
    echo "it under the terms of the GNU General Public License as published by"
    echo "the Free Software Foundation, either version 3 of the License, or"
    echo "(at your option) any later version."
    echo ""
    echo "MoleTools is distributed in the hope that it will be useful,"
    echo "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    echo "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    echo "GNU General Public License for more details."
    echo ""
    echo "You should have received a copy of the GNU General Public License"
    echo "along with MoleTools.  If not, see <http://www.gnu.org/licenses/>."
    echo "*/"
    echo ""
    echo ""
    echo "#include \"Metallo_RF.hpp\""
    echo ""
    echo "std::string Metallo_RF::getTree(unsigned int elem, unsigned int predictor){"
    echo "    //Serialized trees generated from Scikit Learn trees (see printTrees.py and getTrees.bash)"

    # script to extract trees
    echo "    if (predictor == 1){"
    echo "        switch (elem){"
    for ((i=0;i<ntree;i++))
    do
          j=`echo "$i+1" | bc -l`
          echo "      case $i:"
          tree=`awk '{ if(NR=="'"${j}"'") print $0}' data/trees`
          echo "        return \" $tree\";"
          echo "        break;"
    done
    echo "      default:"
    echo "        return \"\";"
    echo "        break;"
	echo "        }"
    echo "    }"
    echo "}"
    echo ""

    echo "unsigned int Metallo_RF::getNTree(){"
    echo "  return $ntree;"
    echo "}"
fi
