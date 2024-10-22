Autocycler's documentation lives in this [repository's wiki](https://github.com/rrwick/Autocycler/wiki). Most users should head there for the latest version of the documentation. This directory is intended for contributors or those who need to build the documentation offline using mdBook.

### Prerequisites:
You must have [mdBook](https://rust-lang.github.io/mdBook/) installed to run the scripts in this directory.

### Running `./build.sh` will:
* Clone the Autocycler wiki
* Create an mdBook book
* Make some necessary modifications to the Markdown files:
  - `create_summary.py` generates a table of contents based on the wiki's `_Sidebar.md`.
  - `fix_links_add_headers.py` fixes internal links and adds necessary headers for the mdBook format.
* Build the mdBook book
* Tarball the result into `autocycler_docs.tar.gz`

After extracting `autocycler_docs.tar.gz`, users can open the documentation by double-clicking `autocycler_docs.html`.
