# optimisation_celery

Demonstration files showing celery being used to manage optimisation tasks.

Start a worker using the command (ideally within a virtualenv)

celery -A optimiser worker --loglevel=INFO -concurrency=<N> -Ofair
