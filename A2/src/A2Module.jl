module A2Module

#Importing A1Module 
include("A1Module.jl") # Module written by me, Kim Laberinto
using .A1Module: SwannsBracketingMethod, GoldenSectionSearch

# Other Imports
using LinearAlgebra #External module for taking norm
using Memento #Invenia Module for logging
using Printf #External module for formatting strings
using ValueHistories #External Package for keeping track of values

end