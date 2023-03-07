
"""
Module for determining point group of a molecule within some tolerance tol.
Includes module CharacterTables which generates information such as symmetry elements,
character tables, multiplication tables, and irreducible representation matrices.
Authors: Stephen Goodlett and Nate Kitzmiller
"""


module Symmetry

using Molecules

tol = 1E-4

include("sea.jl")
include("moit.jl")
include("symmetry_finding.jl")
include("flowchart.jl")
include("CharacterTables/CharacterTables.jl")

end

