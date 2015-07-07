

#' @describeIn Experiment Showing simulated experiment as an image
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_raster scale_fill_gradient2 aes
#' @exportMethod image
#' @param x \link{Experiment-class} object
#' @param ... currently ignored
setMethod(
    "image",
    "Experiment",
    #     signature(object="Experiment"),
    function(x, ...)
    {
        stopifnot(x@simulated)
        int <- x@simFinal
        # log2 transform, zero-centering
        int <- log2(int)
        int <- sweep(int, 1, rowMeans(int, na.rm=T), '-')
        xm <- melt(int)
        xm$peptides <- ordered(xm$peptides, levels=rev(x@proteins@pepNames))
        #
        ggplot(xm, aes(x=samples, y=peptides, fill=value)) +
            geom_raster() +
            scale_fill_gradient2()
    }
)


