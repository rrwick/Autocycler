# Clean any previous build:
rm -rf Autocycler.wiki src book.toml book autocycler_docs.tar.gz autocycler_docs.html autocycler_docs/

# Clone the repo:
git clone https://github.com/rrwick/Autocycler.wiki.git

# Create an mdbook:
mdbook init --force

# Set some build options:
echo "\n[output.html]" >> book.toml
echo "no-section-label = true" >> book.toml
echo "preferred-dark-theme = \"ayu\"" >> book.toml

# Move pages to mdbook src:
cp Autocycler.wiki/*.md src/
rm src/_Sidebar.md

# Change image locations to relative paths:
for f in src/*.md; do
    sed -i 's|https://github.com/rrwick/Autocycler/wiki/||g' "$f"
done

# Fix inter-page links and add headers:
for f in src/*.md; do
    ./fix_links_add_headers.py "$f"
done

# Convert the sidebar page into a summary page:
./create_summary.py Autocycler.wiki/_Sidebar.md > src/SUMMARY.md

# Build the book:
mdbook build

# Copy images into the book directory:
cp -r Autocycler.wiki/images book

# Add a bit more space above headers in the TOC:
sed -i 's/margin: 5px 0px;/margin: 20px 0px 5px 0px;/' book/css/general.css

# Set up for easier distribution:
mv book autocycler_docs
mv autocycler_docs/index.html autocycler_docs.html
sed -i 's|href="|href="autocycler_docs/|g' autocycler_docs.html
sed -i 's|src="|src="autocycler_docs/|g' autocycler_docs.html
sed -i 's|var path_to_root = "";|var path_to_root = "autocycler_docs/";|' autocycler_docs.html
tar -zcvf autocycler_docs.tar.gz --owner=0 --group=0 autocycler_docs.html autocycler_docs/

# Clean up:
rm -rf Autocycler.wiki src book.toml __pycache__
