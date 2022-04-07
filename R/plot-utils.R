.add_category_label <- function(data, x_cat, y_cat, cutoff, categories){
    data$Category <- categories[5L]
    x_only <- abs(data[,x_cat]) >= cutoff &
        abs(data[,y_cat]) < cutoff
    data$Category[x_only] <- categories[1L]
    y_only <- abs(data[,y_cat]) >= cutoff &
        abs(data[,x_cat]) < cutoff
    data$Category[y_only] <- categories[2L]
    x_and_y <- abs(data[,x_cat]) >= cutoff &
        abs(data[,y_cat]) >= cutoff
    data$Category[x_and_y] <- categories[3L]
    x_opposite_y <- abs(data[,x_cat]) >= cutoff &
        abs(data[,y_cat]) >= cutoff &
        sign(data[,x_cat]) != sign(data[,y_cat])
    data$Category[x_opposite_y] <- categories[4L]
    data
}

.plot_scatter <- function(data, x, y, categories, cutoff, xlim, ylim){
    data <- .add_category_label(data,
                                x_cat = x,
                                y_cat = y,
                                cutoff = cutoff,
                                categories = names(categories))
    plot_out <- ggplot(as.data.frame(data),
                       aes_string(x = x, y = y)) +
        geom_hline(yintercept = cutoff, colour = "gray") +
        geom_hline(yintercept = -cutoff, colour = "gray") +
        geom_vline(xintercept = cutoff, colour = "gray") +
        geom_vline(xintercept = -cutoff, colour = "gray") +
        geom_point(aes_string(colour = "Category"), shape = 19L, alpha = 0.65) +
        scale_x_continuous(limits = xlim) +
        scale_y_continuous(limits = ylim) +
        scale_colour_manual(values = categories) +
        theme_bw()
    plot_out
}
