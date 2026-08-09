#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Train an XGBoost model on the HCC training set and predict mutations in a new sample set.

This is a cleaned-up, parameterized version of 预测集验证.py. The output only
contains the mutation ID and the prediction results (predicted class and probability).

Usage:
    python3 predict.py \
        --train ../data/HCC_training_all.txt \
        --predict feature_add10.txt \
        --output pred-results.txt \
        [--sample-name Figure7]
"""

import argparse
import os

import numpy as np
import pandas as pd
import xgboost as xgb


MODEL_PARAMS = {
    "eta": 0.01,
    "max_depth": 6,
    "objective": "binary:logistic",
    "eval_metric": "logloss",
    "gamma": 1.5,
    "subsample": 0.9,
    "colsample_bytree": 0.25,
    "n_estimators": 189,
    "min_child_weight": 2.7,
    "lambda": 0.5,
    "alpha": 0.5,
    "random_state": 68,
}


def clean_column_names(df):
    """Strip whitespace from column names and convert them to lowercase."""
    df.columns = df.columns.str.strip()
    df.columns = df.columns.str.replace(r"\s+", "", regex=True)
    df.columns = df.columns.str.lower()
    return df


def read_data(file_path, data_name):
    """Read a tab-separated table and check the required columns."""
    try:
        df = pd.read_csv(file_path, delimiter="\t", encoding="utf-8")
    except Exception:
        try:
            df = pd.read_csv(file_path, delimiter=r"\s+", encoding="gbk")
        except Exception:
            df = pd.read_csv(file_path, delimiter=r"\s+", encoding="utf-8-sig")

    df = clean_column_names(df)

    if "train" in data_name.lower():
        required_cols = ["context", "type", "region", "label"]
    else:
        required_cols = ["context", "type", "region"]

    actual_cols = df.columns.tolist()
    missing_cols = [col for col in required_cols if col not in actual_cols]
    if missing_cols:
        print("\n%s actual columns (first 20):" % data_name)
        for i, col in enumerate(actual_cols[:20]):
            print("  %d. '%s'" % (i + 1, col))
        raise ValueError("%s is missing required columns: %s" % (data_name, missing_cols))

    print("%s read successfully: shape=%s, required columns present" % (data_name, df.shape))
    return df


def convert_to_dummies(data, cat_dict=None, is_train=True):
    """One-hot encode context/type/region. Prediction data uses the training
    dummy columns so that train and test feature sets are identical."""
    dummy_cols = []
    for col in ["context", "type", "region"]:
        if is_train:
            dummies = pd.get_dummies(data[col], prefix=col, dummy_na=True)
            dummy_cols.extend(dummies.columns.tolist())
        else:
            if col not in data.columns:
                print("Prediction data is missing '%s', creating empty column" % col)
                data[col] = np.nan
            train_dummy_cols = [c for c in cat_dict.keys() if c.startswith(col + "_")]
            dummies = pd.get_dummies(data[col], prefix=col, dummy_na=True)
            for train_col in train_dummy_cols:
                if train_col not in dummies.columns:
                    dummies[train_col] = 0
            dummies = dummies[train_dummy_cols]

        data = pd.concat([data, dummies], axis=1)
        if col in data.columns:
            data = data.drop(col, axis=1)
    return data, dummy_cols if is_train else None


def preprocess_data(df, cat_columns_dict=None, is_train=True):
    if is_train:
        df, dummy_cols = convert_to_dummies(df)
        cat_columns_dict = {col: None for col in dummy_cols}
        drop_cols = [c for c in ["id", "label", "conflict_num", "mutation_frequency"] if c in df.columns]
        X = df.drop(drop_cols, axis=1)
        y = df["label"]
        print("Training set after preprocessing: X=%s, y=%s" % (X.shape, y.shape))
        return X, y, cat_columns_dict
    else:
        df, _ = convert_to_dummies(df, cat_dict=cat_columns_dict, is_train=False)
        drop_cols = [c for c in ["id", "conflict_num", "mutation_frequency"] if c in df.columns]
        X = df.drop(drop_cols, axis=1)
        print("Prediction set after preprocessing: X=%s" % (X.shape,))
        return X, df


def main():
    parser = argparse.ArgumentParser(description="Train XGBoost on HCC training data and predict new mutations.")
    parser.add_argument("--train", required=True, help="training set (e.g. HCC_training_all.txt)")
    parser.add_argument("--predict", required=True, help="prediction feature file (e.g. feature_add10.txt)")
    parser.add_argument("--output", required=True, help="output prediction result file")
    parser.add_argument("--sample-name", default="prediction", help="sample/cohort name (used in the log only)")
    args = parser.parse_args()

    print("=" * 60)
    print("Step 1: reading training set")
    train_df = read_data(args.train, "训练集")
    X_train, y_train, cat_columns = preprocess_data(train_df, is_train=True)

    print("=" * 60)
    print("Step 2: reading prediction set")
    pred_df = read_data(args.predict, "预测集")
    X_pred, raw_pred_df = preprocess_data(pred_df, cat_columns_dict=cat_columns, is_train=False)

    # align prediction features with training features
    missing_features = set(X_train.columns) - set(X_pred.columns)
    for feat in missing_features:
        X_pred[feat] = 0
    X_pred = X_pred[X_train.columns]

    print("Final feature columns: training=%d, prediction=%d" % (len(X_train.columns), len(X_pred.columns)))
    if list(X_train.columns) != list(X_pred.columns):
        raise ValueError("Feature columns between training and prediction sets are inconsistent!")
    print("Feature columns are identical. Training model...")

    print("=" * 60)
    print("Step 3: training XGBoost model (%s)" % args.sample_name)
    model = xgb.XGBClassifier(**MODEL_PARAMS)
    model.fit(X_train, y_train)

    print("Step 4: predicting mutations")
    pred_label = model.predict(X_pred)
    pred_prob = model.predict_proba(X_pred)[:, 1]

    # output: mutation ID + prediction result only
    result = pd.DataFrame({
        "id": raw_pred_df["id"],
        "pred_label": pred_label,
        "pred_prob": pred_prob,
    })

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, encoding="utf-8")

    n_total = len(result)
    n_mut = int((result["pred_label"] == 1).sum())
    n_wt = n_total - n_mut
    print("\nPrediction finished: %d mutations, %d predicted as mutated (1), %d predicted as wild-type (0)" %
          (n_total, n_mut, n_wt))
    print("Prediction results saved to: %s" % args.output)


if __name__ == "__main__":
    main()
