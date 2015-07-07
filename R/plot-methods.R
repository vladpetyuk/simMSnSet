
#' @export
setGeneric("plot", function(x, y, z, ...) standardGeneric("plot"))


#' @importFrom ggplot2 geom_line geom_point
setMethod(
    "plot",
    signature("Experiment", "character", "ANY"),
    definition=function(x, y, z, ...){
        idx <- x@proteins@features$prot == y
        m <- log10(t(x@simFinal[idx,]))
        mm <- melt(m)
        p <- ggplot(mm, aes(x=samples, y=value, group=peptides, color=peptides)) +
            geom_line(size=1.5) +
            geom_point(size=2, shape=21, fill='white')
        suppressWarnings(graphics::plot(p))
    }
)


#' @importMethodsFrom Biobase fData exprs
setMethod(
    "plot",
    signature("MSnSet", "character", "ANY"),
    definition=function(x, y, z, protColName="prot", relative=FALSE, ...){
        idx <- fData(x)[[protColName]] == y
        m <- log10(exprs(x)[idx,])
        if(relative)
            m <- sweep(m, 1, rowMeans(m, na.rm=TRUE), '-')
        mm <- melt(m)
        p <- ggplot(mm, aes(x=samples, y=value, group=peptides, color=peptides)) +
            geom_line(size=1.5) +
            geom_point(size=2, shape=21, fill='white')
        suppressWarnings(graphics::plot(p))
    }
)




setMethod(
    "plot",
    signature("MSnSet", "MSnSet", "character"),
    definition=function(x, y, z, protColName="prot", ...){
        idx <- fData(x)[[protColName]] == z
        mpep <- log10(exprs(x)[idx,])
        idx <- fData(y)[[protColName]] == z
        mpro <- log10(exprs(y)[idx,])
        mpep <- sweep(mpep, 1, rowMeans(mpep, na.rm=TRUE), '-')
        mpro <- sweep(mpro, 1, rowMeans(mpro, na.rm=TRUE), '-')
        mmpep <- melt(mpep)
        mmpro <- melt(mpro)
#         p <- ggplot(mm, aes(x=samples, y=value, group=peptides, color=peptides)) +
#             geom_line(size=1.5) +
#             geom_point(size=2, shape=21, fill='white')
#         suppressWarnings(plot(p))
    }
)
