#Openness to Experience----
#' Four Inventories of Openness to Experience
#'
#' A response matrix (\emph{n} = 794) of all four Openness to Experience
#' inventories from Christensen, Cotter, & Silvia (2019). The
#' key provides inventory, facet, and item description information
#' for the item labels. Note that because of NEO's copyrights the
#' items have been shortened and paraphrased
#'
#' @name openness
#'
#' @aliases openness
#' @aliases openness.key
#'
#' @docType data
#'
#' @usage data(openness)
#' @usage data(openness.key)
#'
#' @format A 794 x 138 response matrix (openness) and 138 x 7 matrix (openness.key).
#' Here are detailed descriptions of the key:
#'
#' \itemize{
#'
#' \item \code{Inventory} --- The personality inventory the item belongs to
#'
#' \item \code{Facet} --- The personality inventory defined facet
#'
#' \item \code{JPA.Domains} --- The broad domains identified by Christensen, Cotter, and Silvia (2019)
#'
#' \item \code{JPA.Facets} --- The facets identified by Christensen, Cotter, and Silvia (2019)
#'
#' \item \code{Item.Label} --- The labels used in Christensen, Cotter, and Silvia (2019)
#'
#' \item \code{Item.Description} --- Descriptions of each item. Note that the NEO-PI-3 items are protected by
#' copyright and therefore have been paraphrased. These item descriptions
#' do not represent the item as given to the participant
#'
#' \item \code{Reversed} --- Whether an item should be reversed or not (\code{openness} is already reversed)
#'
#' }
#'
#'
#' @keywords datasets
#'
#' @references
#' Christensen, A. P., Cotter, K. N., & Silvia, P. J. (2019).
#' Reopening openness to experience: A network analysis of four openness to experience inventories.
#' \emph{Journal of Personality Assessment}, \emph{101}, 574-588.
#'
#' @examples
#' # Loading data
#' data("openness")
#' data("openness.key")
#'
#' # Change item labels
#' colnames(openness) <- openness.key$Item.Description
#'
NULL
#----