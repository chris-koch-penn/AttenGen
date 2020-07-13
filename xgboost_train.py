# Chris Koch, 2020. Model adapted from VaxignML by Edison Ong.
import numpy as np
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import GridSearchCV, StratifiedKFold
import joblib
from pathlib import Path
import pickle
import sys
from sklearn.preprocessing import MinMaxScaler


def load_victors(base):
    uniprot = pickle.load(open(base / "victors_uniprot.pkl", "rb"))
    vf_bacteria = pickle.load(open(base / "vf_bacteria.pkl", "rb"))
    vf_viruses = pickle.load(open(base / "vf_viruses.pkl", "rb"))
    return uniprot, vf_bacteria, vf_viruses


def load_protegen(base):
    uniprot = pickle.load(open(base / "protegen_uniprot.pkl", "rb"))
    vf_bacteria = pickle.load(open(base / "protegen_bacteria.pkl", "rb"))
    vf_viruses = pickle.load(open(base / "protegen_viruses.pkl", "rb"))
    return uniprot, vf_bacteria, vf_viruses


def train_model(IS_VICTORS):
    # Load data.
    print("Loading data...")
    base = Path("./data")
    uniprot, vf_bacteria, vf_viruses = load_protegen(
        base) if IS_VICTORS else load_victors(base)
    x_train, y_train = [], []
    x_train = uniprot["x_train"] + \
        vf_bacteria["x_train"] + vf_viruses["x_train"]
    y_train = uniprot["y_train"] + \
        vf_bacteria["y_train"] + vf_viruses["y_train"]
    x_train, y_train = np.asarray(x_train), np.asarray(y_train)

    # Declare model and pipeline.
    train_on = sys.argv[1] if len(sys.argv) > 1 else "gpu"
    est = None
    if train_on == "cpu":
        est = XGBClassifier(objective='binary:logistic', nthread=15,
                            eval_metric='auc', random_state=26)
    else:
        est = XGBClassifier(objective='binary:logistic', nthread=15,
                            eval_metric='auc', random_state=26,
                            tree_method='gpu_hist', predictor='gpu_predictor')
    estPipe = Pipeline(
        [('feature_selection', SelectKBest()), ('classification', est)])
    grid = [{
        # "feature_selection__k": [110, 130, 150, 180, 210, 240],
        "feature_selection__k": [240, 280, 320, 360, 400, 440],
        'classification__learning_rate': [0.3, 0.1],
        'classification__n_estimators': [80, 100, 120, 140, 160, 300, 1000],
        'classification__max_depth': [3, 6, 9],
        'classification__min_child_weight': [1, 3],
        'classification__scale_pos_weight': [1, 6],
        'classification__max_delta_step': [0, 3],
    }]

    # Train model.
    print("Training model...")
    prefix = "victors" if IS_VICTORS else "protegen"
    cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=12)
    xgb = GridSearchCV(estimator=estPipe, param_grid=grid,
                       cv=cv, verbose=1)
    xgb.fit(x_train, y_train)
    y_prob = xgb.predict_proba(x_train)
    joblib.dump(xgb, f"./saved_models/{prefix}_xgboost_model.joblib")
    joblib.dump(y_prob, f"./saved_models/{prefix}_xgboost_scores.joblib")


if __name__ == "__main__":
    IS_VICTORS, IS_PROTEGEN = True, False
    train_model(IS_VICTORS)
    # train_model(IS_PROTEGEN)
