if (getRversion() >= "2.15.1")
    utils::globalVariables(c("sequence", "gonnet"))

# set const number of amino acids plus gap
AA_COUNT = length(append(Biostrings::AA_STANDARD, "-"))