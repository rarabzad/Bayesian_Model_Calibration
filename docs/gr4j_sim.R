gr4j_sim <- function(P_xts, 
                     T_xts,
                     PET_xts = NULL,
                     params,
                     area_km2,
                     lat = NULL,
                     return_states = FALSE)
{
  # Load libraries
  if (!requireNamespace("xts", quietly = TRUE)) stop("xts package required")
  if (!requireNamespace("zoo", quietly = TRUE)) stop("zoo package required")
  
  # Align and merge input time series
  if (is.null(PET_xts)) {
    data_all <- merge(P_xts, T_xts)
  } else {
    data_all <- merge(P_xts, T_xts, PET_xts)
  }
  data_all <- na.omit(data_all)
  
  # Extract input vectors
  P <- as.numeric(data_all[, 1])
  T <- as.numeric(data_all[, 2])
  dates <- index(data_all)
  n <- length(dates)
  
  # Handle PET (compute if not given)
  if (is.null(PET_xts)) {
    if (is.null(lat)) stop("Latitude must be provided when PET_xts is NULL.")
    JD <- as.numeric(strftime(dates, "%j"))
    if (requireNamespace("airGR", quietly = TRUE)) {
      PET <- airGR::PE_Oudin(JD = JD, Temp = T, Lat = lat, LatUnit = "deg")
    } else stop("airGR package required for Oudin PET calculation.")
  } else {
    PET <- as.numeric(data_all[, 3])
  }
  
  # Unpack parameters
  X1 <- params[1]; X2 <- params[2]; X3 <- params[3]; X4 <- params[4]
  TT <- params[5]; DDF <- params[6]
  
  # Initialize states
  S_prev <- 0.5 * X1
  SWE_prev <- 0
  S <- numeric(n)
  SWE <- numeric(n)
  U_eff <- numeric(n)
  
  for (t in seq_len(n)) {
    if (T[t] < TT) {
      SWE_cur <- SWE_prev + P[t]
      Rain_in <- 0
    } else {
      melt_pot <- DDF * (T[t] - TT)
      melt_act <- min(melt_pot, SWE_prev)
      SWE_cur <- SWE_prev - melt_act
      Rain_in <- P[t] + melt_act
    }
    
    # Net rainfall and evapotranspiration
    Pn <- max(Rain_in - PET[t], 0)
    En <- max(PET[t] - Rain_in, 0)
    St_ratio <- S_prev / X1
    
    if (Pn > 0) {
      Ps <- X1 * (1 - St_ratio^2) * tanh(Pn / X1) / (1 + St_ratio * tanh(Pn / X1))
      Es <- 0
    } else {
      Ps <- 0
      Es <- S_prev * (2 - St_ratio) * tanh(En / X1) / (1 + (1 - St_ratio) * tanh(En / X1))
    }
    
    S_cur <- S_prev - Es + Ps
    perc <- S_cur * (1 - (1 + ((4/9)*(S_cur/X1))^4)^(-0.25))
    S_new <- S_cur - perc
    U_eff[t] <- perc + (Pn - Ps)
    
    # Update states
    S_prev <- S_new
    SWE_prev <- SWE_cur
    S[t] <- S_new
    SWE[t] <- SWE_cur
  }
  
  # Unit hydrographs
  n1 <- ceiling(X4); n2 <- ceiling(2 * X4)
  t1 <- seq(1, n1); SH1 <- pmin((t1 / X4)^2.5, 1)
  SH2 <- numeric(n2)
  for (i in 1:n2) {
    if (i <= X4) {
      SH2[i] <- 0.5 * (i / X4)^2.5
    } else if (i <= 2 * X4) {
      SH2[i] <- 1 - 0.5 * ((2 - i / X4)^2.5)
    } else {
      SH2[i] <- 1
    }
  }
  UH1 <- diff(c(0, SH1))
  UH2 <- diff(c(0, SH2))
  
  # Routing: convolution
  Q9 <- stats::filter(c(rep(0, length(UH1)), 0.9 * U_eff), filter = UH1, sides = 1)[-(1:length(UH1))]
  Q1 <- stats::filter(c(rep(0, length(UH2)), 0.1 * U_eff), filter = UH2, sides = 1)[-(1:length(UH2))]
  
  # Routing reservoir
  Qr <- Qd <- R_res <- numeric(n)
  R_prev <- 0
  for (t in seq_len(n)) {
    F <- X2 * (R_prev / X3)^3.5
    R_temp <- max(0, R_prev + Q9[t] + F)
    Qr[t] <- R_temp * (1 - (1 + (R_temp / X3)^4)^(-0.25))
    R_res[t] <- R_temp - Qr[t]
    Qd[t] <- max(0, Q1[t] + F)
    R_prev <- R_res[t]
  }
  
  # Convert to discharge in m³/s
  Qsim_mm <- Qr + Qd
  Qsim_m3s <- Qsim_mm * area_km2 * 1e3 / 86400
  Qsim_xts <- xts::xts(Qsim_m3s, order.by = dates)
  
  if (!return_states) {
    return(Qsim_xts)
  } else {
    SWE_xts <- xts::xts(SWE, order.by = dates)
    S_xts <- xts::xts(S, order.by = dates)
    PET_xts_out <- xts::xts(PET, order.by = dates)
    return(list(Q = Qsim_xts, SWE = SWE_xts, S = S_xts, PET = PET_xts_out))
  }
}
