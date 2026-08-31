# pandoc "Understanding the Gillespie Algorithm.md" \
#   -o "Understanding the Gillespie Algorithm.html" \
#   --standalone \
#   --mathjax \
#   --embed-resources \
#   --from=markdown+markdown_in_html_blocks+tex_math_dollars-implicit_figures

# pandoc "Understanding the Gillespie Algorithm.md" \
#   -o "Understanding the Gillespie Algorithm.html" \
#   --standalone \
#   --mathml

pandoc "Understanding the Gillespie Algorithm.md" \
  -o "Understanding the Gillespie Algorithm.html" \
  --standalone \
  --mathjax=https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg-full.js \
  --embed-resources \
  --from=markdown+markdown_in_html_blocks+tex_math_dollars-implicit_figures