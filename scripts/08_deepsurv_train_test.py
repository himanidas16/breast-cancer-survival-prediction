#!/usr/bin/env python3
"""
08_deepsurv_train_test.py
Clean DeepSurv train/test benchmark.

Expected local/private files:
  data/deepsurv_export_LOCAL_ONLY/deepsurv_expr_LOCAL_ONLY.csv
  data/deepsurv_export_LOCAL_ONLY/deepsurv_clinical_LOCAL_ONLY.csv

Do not upload the patient-level CSVs or prediction outputs.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import torch
import torchtuples as tt
from pycox.models import CoxPH
from pycox.evaluation import EvalSurv

ROOT = Path(".")
EXPORT_DIR = ROOT / "data" / "deepsurv_export_LOCAL_ONLY"
TABLE_DIR = ROOT / "results" / "summary_tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

EXPR_CSV = EXPORT_DIR / "deepsurv_expr_LOCAL_ONLY.csv"
CLIN_CSV = EXPORT_DIR / "deepsurv_clinical_LOCAL_ONLY.csv"

if not EXPR_CSV.exists() or not CLIN_CSV.exists():
    raise FileNotFoundError("Run scripts/07_export_for_deepsurv.R first and keep exported CSVs local/private.")

expr = pd.read_csv(EXPR_CSV)
clin = pd.read_csv(CLIN_CSV)

if not expr["sample_id"].equals(clin["sample_id"]):
    raise ValueError("Sample order mismatch between expression and clinical files.")

X = expr.drop(columns=["sample_id"]).values.astype("float32")
time = clin["os_time"].values.astype("float32")
event = clin["os_event"].values.astype("int64")
train_mask = clin["set"].values == "train"
test_mask = clin["set"].values == "test"

scaler = StandardScaler()
X_train = scaler.fit_transform(X[train_mask]).astype("float32")
X_test = scaler.transform(X[test_mask]).astype("float32")

y_train = (time[train_mask], event[train_mask])
y_test = (time[test_mask], event[test_mask])

in_features = X_train.shape[1]
net = tt.practical.MLPVanilla(
    in_features=in_features,
    num_nodes=[128, 64],
    out_features=1,
    batch_norm=True,
    dropout=0.2,
    activation=torch.nn.ReLU,
)

model = CoxPH(net, tt.optim.Adam(lr=1e-3))
log = model.fit(
    X_train,
    y_train,
    batch_size=256,
    epochs=100,
    callbacks=[tt.callbacks.EarlyStopping(patience=10)],
    verbose=True,
    val_data=(X_test, y_test),
)

model.compute_baseline_hazards()
# Pycox concordance_td is not exactly Harrell C-index, but it is a practical DeepSurv eval metric.
surv = model.predict_surv_df(X_test)
ev = EvalSurv(surv, time[test_mask], event[test_mask], censor_surv="km")
c_td = ev.concordance_td("antolini")

summary = pd.DataFrame([
    {
        "model": "DeepSurv",
        "train_patients": int(train_mask.sum()),
        "test_patients": int(test_mask.sum()),
        "train_events": int(event[train_mask].sum()),
        "test_events": int(event[test_mask].sum()),
        "input_genes": int(in_features),
        "hidden_layers": "128,64",
        "dropout": 0.2,
        "batch_size": 256,
        "max_epochs": 100,
        "cindex_td_antolini": float(c_td),
    }
])
summary.to_csv(TABLE_DIR / "deepsurv_summary_generated.csv", index=False)

# Save model locally only.
LOCAL_MODEL_DIR = ROOT / "data" / "processed"
LOCAL_MODEL_DIR.mkdir(parents=True, exist_ok=True)
torch.save(model.net.state_dict(), LOCAL_MODEL_DIR / "deepsurv_weights_LOCAL_ONLY.pt")

print("DeepSurv complete. Aggregate summary saved. Patient-level outputs were not written.")
