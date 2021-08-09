plot_km <- function(data,
                    factor = "1",
                    factor_labels = NULL,
                    time_breaks = 1,
                    confidence_interval = TRUE,
                    colors = NULL,
                    risk_table = TRUE,
                    time_range = c(0, 20),
                    title = "",
                    ggtheme = ggplot2::theme_light(),
                    table_height = 0.35,
                    ylim = c(0, 1),
                    ofn = NULL,
                    width = 8,
                    height = 10,
                    ggtheme_table = survminer::theme_cleantable(base_size = 6, base_family = "Arial"),
                    ...) {
    data <- as_tibble(data)

    fit <- survminer::surv_fit(as.formula(paste0("survival::Surv(y, death)~", factor)), data = data)

    if (is.null(factor_labels) && factor == "1") {
        factor_labels <- ""
    }

    data <- as.data.frame(data)

    g <- survminer::ggsurvplot(fit,
        data = data,
        break.time.by = time_breaks,
        conf.int = confidence_interval,
        palette = colors,
        title = title,
        risk.table = risk_table,
        xlim = time_range,
        ggtheme = ggtheme,
        legend.labs = factor_labels,
        tables.height = table_height,
        tables.theme = ggtheme_table,
        ylim = ylim,
        ...
    )

    if (!is.null(ofn)) {
        ggsave(file = ofn, print(g), width = width, height = height)
    }

    return(g)
}


