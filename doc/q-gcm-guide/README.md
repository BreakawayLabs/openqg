## Building

To build the documentation run the following command:

    pdflatex q-gcm.tex; pdflatex q-gcm.tex

This will generate a file called q-gcm.pdf in the current directory. You need to run pdflatex twice so as to correctly generate the table of contents.

## Dependencies

To build the documentation you must have `pdflatex` installed. On Ubuntu systems this can be installed with the following command:

    sudo apt-get install texlive-latex-extra
