#! /bin/bash
echo
echo "Fixing LE tables"
sed -i -z 's/\n)/)/g' data/GSEA/LE_*
sed -i -z 's/\n,/,/g' data/GSEA/LE_*
sed -i -z 's/\n"/"/g' data/GSEA/LE_*
echo
echo "Tables were fixed!"
echo

