using Documenter
using PriceIndexes

makedocs(
    sitename = "PriceIndexes",
    doctest = true,
    format = Documenter.HTML(sidebar_sitename = false),
    modules = [PriceIndexes],
    checkdocs=:exports #add to try and fix error
)



# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
