import joblib
import pickle
from xgboost_train import load_victors
import numpy as np
from pathlib import Path
import random
from sklearn.metrics import classification_report, confusion_matrix


def load_data(p1):
    # Load data and models.
    xgb = joblib.load(p1)
    print("Best parameters: ", xgb.best_params_)
    print("Mean Weighted F1 score on 3-fold cross-validation for best estimator: ", xgb.best_score_)


if __name__ == "__main__":
    protegen_p1 = "./saved_models/protegen_xgboost_model.joblib"
    vf_p1 = "./saved_models/victors_xgboost_model.joblib"
    print("VIRAL MODEL: ")
    load_data(vf_p1)
    print("PROTEGEN MODEL: ")
    load_data(protegen_p1)
