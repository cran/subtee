get_prca_data <- function() {
  # data_url = "https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/pdf/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/adv_prostate_ca.zip"
  data_url = "https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/adv_prostate_ca.zip"
  if (httr::http_error(data_url)) { # network is down = message (not an error)
    message("No internet connection or data source broken.")
    return(NULL)
  } else { # network is up = proceed to download via curl
    message("Downloading remote dataset.")
    temp <- tempfile()
    utils::download.file(data_url, temp, cacheOK = FALSE)
    prca = utils::read.csv(unz(temp, "adv_prostate_ca/adv_prostate_ca.csv"))
    names(prca) = toupper(names(prca))
    return(prca)
  }
}
