library(future)
library(progressr)

future::plan(future::multisession)
progressr::handlers(global = TRUE)
progressr::handlers(
  progressr::handler_txtprogressbar(),
  progressr::handler_progress(format = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta")
)
