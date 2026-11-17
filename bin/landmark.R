# ============================================================================
# landmark.R: survival times counted from the end of induction (a landmark), not from diagnosis
# ============================================================================
# Sourced by tractor_genesis_adapter.R; tested in tests/test_tractor_model.R.
#
# Decided 17 Nov 2026: EFS and OS start at the end of induction, because the survival models are
# stratified on (or adjusted for) end-of-induction MRD, which does not exist before then. Counting
# from diagnosis would credit every child with the event-free weeks of induction it had to survive
# to have an MRD value at all (immortal time).
#
# The handoff bundle ships each time in days from diagnosis (<trait>_time) and each child's
# end-of-induction day, also in days from diagnosis (eoi_day: day 29 in COG and REDIAL, day 46 at
# St. Jude, or the child's own date where it is known). The landmark time is time - eoi_day. A child
# whose event or censoring falls AT OR BEFORE eoi_day was not at risk at the landmark, and a child
# without eoi_day cannot be placed on the new clock: both are left out of that trait and counted.
#
# landmark_eoi(time, event, eoi)
#   time, event: the trait's columns as shipped (days from diagnosis; 0/1); eoi: eoi_day
#   returns list(time      = time - eoi, NA where the child is left out (or time / event is NA),
#                excluded  = TRUE for a child left out by the landmark (time and event known),
#                n_before  = how many of them had their event or censoring at or before eoi_day,
#                n_no_eoi  = how many had no eoi_day)
# Rows whose time or event is NA are left alone here: the adapter's complete-case filter drops and
# counts them, as for every trait.
# ============================================================================
landmark_eoi <- function(time, event, eoi) {
    time <- as.numeric(time); eoi <- suppressWarnings(as.numeric(eoi))
    known <- !is.na(time) & !is.na(event)
    no_eoi <- known & is.na(eoi)
    before <- known & !is.na(eoi) & time <= eoi
    excluded <- no_eoi | before
    out <- time - eoi
    out[excluded | !known] <- NA_real_
    list(time = out, excluded = excluded, n_before = sum(before), n_no_eoi = sum(no_eoi))
}

# parse_sample_subset("study_group=REDIAL") -> list(column = "study_group", values = "REDIAL");
# several values: "study_group=REDIAL,SJ_TXV". Stops, naming the problem, on anything else.
parse_sample_subset <- function(spec) {
    m <- regmatches(spec, regexec("^([A-Za-z][A-Za-z0-9_.]*)=(.+)$", spec))[[1]]
    if (length(m) != 3) stop("--sample_subset '", spec, "': expected column=value[,value...] (e.g. study_group=REDIAL)")
    vals <- trimws(strsplit(m[3], ",")[[1]])
    vals <- vals[vals != ""]
    if (length(vals) == 0) stop("--sample_subset '", spec, "': no value after '='")
    list(column = m[2], values = vals)
}
