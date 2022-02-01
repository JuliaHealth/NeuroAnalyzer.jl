#' @title jaccard_similarity
#' 
#' @description Calculates Jaccard similarity between two vectors
#' 
#' @param x The vector to analyze
#' @param y The vector to analyze
#' @export
#' 
#' @return

jaccard_similarity <- function(x, y)
{
    intersection = length(intersect(x, y))
    union = length(x) + length(y) - intersection
    intersection / union
}