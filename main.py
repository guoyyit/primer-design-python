"""Entry point: uvicorn main:app --host 127.0.0.1 --port 8000"""
import os

from primer_design.api import create_app

app = create_app(os.environ.get("PRIMER_DESIGN_HOME", ".primer_design_data"))
