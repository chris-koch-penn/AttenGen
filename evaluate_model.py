import joblib
import pickle
from xgboost_train import load_victors
from genetic_algorithm import run_on_cpu
import numpy as np
from pathlib import Path
import random
from sklearn.metrics import classification_report, confusion_matrix


def write_results(file_name, y_test, y_pred_test):
    with open(file_name, "w") as metrics:
        metrics.write(classification_report(y_test, y_pred_test))
        c = confusion_matrix(y_test, y_pred_test)
        metrics.write("\nConfusion Matrix (Non-virulent, Virulent): \n")
        metrics.write(np.array2string(c, separator=', '))


def load_data(p1, p2):
    # Load data and models.
    base = Path("./data")
    uniprot, bacteria, viruses = load_victors(base)
    xgb, scores = joblib.load(p1), joblib.load(p2)
    print("Best parameters: ", xgb.best_params_)
    print("Best AUC on 5 fold cross-validation: ", xgb.best_score_)
    return xgb, uniprot, bacteria, viruses


def test_model(xgb, uniprot, bacteria, viruses, name):
    # Evaluate model.
    x_test_uniprot = np.asarray(uniprot["x_test"])
    x_test_bacteria = np.asarray(bacteria["x_test"])
    x_test_viruses = np.asarray(viruses["x_test"])
    v_len = len(x_test_viruses)
    b_len = len(x_test_bacteria)
    balanced_uniprot1 = random.sample(list(x_test_uniprot), v_len)
    balanced_uniprot2 = random.sample(list(x_test_uniprot), b_len)
    x_uni = np.concatenate((balanced_uniprot1, x_test_viruses), axis=0)
    x_uni2 = np.concatenate((balanced_uniprot2, x_test_bacteria), axis=0)
    y_uni = np.concatenate(
        (uniprot["y_test"][:b_len], bacteria["y_test"]), axis=0)
    y_uni2 = np.concatenate(
        (uniprot["y_test"][:v_len], viruses["y_test"]), axis=0)
    y_pred_viruses = xgb.predict(x_uni)
    y_pred_bacteria = xgb.predict(x_uni2)

    # Write results
    write_results(f"xgboost_{name}_bacteria.txt", y_uni, y_pred_bacteria)
    write_results(f"xgboost_{name}_viruses.txt", y_uni2, y_pred_viruses)


if __name__ == "__main__":
    protegen_p1 = "./saved_models/protegen_xgboost_model.joblib"
    protegen_p2 = "./saved_models/protegen_xgboost_model.joblib"
    vf_p1 = "./saved_models/victors_xgboost_model.joblib"
    vf_p2 = "./saved_models/victors_xgboost_scores.joblib"
    xgb, uniprot, bacteria, viruses = load_data(vf_p1, vf_p2)
    test_model(xgb, uniprot, bacteria, viruses, "vf")
    # xgb, uniprot, bacteria, viruses = load_data(protegen_p1, protegen_p2)
    # test_model(xgb, uniprot, bacteria, viruses, "protegen")

# import joblib
# import pickle
# from xgboost_train import load_victors
# import numpy as np
# from pathlib import Path
# import random
# from sklearn.metrics import classification_report, confusion_matrix


# def load_data(p1):
#     # Load data and models.
#     xgb = joblib.load(p1)
#     print("Best parameters: ", xgb.best_params_)
#     print("Mean Weighted F1 score on 3-fold cross-validation for best estimator: ", xgb.best_score_)


# if __name__ == "__main__":
#     protegen_p1 = "./saved_models/protegen_xgboost_model.joblib"
#     vf_p1 = "./saved_models/victors_xgboost_model.joblib"
#     print("VIRAL MODEL: ")
#     load_data(vf_p1)
#     print("PROTEGEN MODEL: ")
#     load_data(protegen_p1)
