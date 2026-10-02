list(years = 1986:2018, ages = 2:11, N_settings = list(process = "off",
    init = "exp"), F_settings = list(process = "rw", mu_form = NULL),
    M_settings = list(process = "rw", mu_form = NULL, mu_supplied = ~M_assumption,
        age_breaks = c(2, 5, 9, 11), first_dev_year = 1986L),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 +
        survey, fill_missing = FALSE))
