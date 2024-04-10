import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, RepeatedKFold, GridSearchCV
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import matplotlib.pyplot as plt

class iLGSScorePredictor:
    def __init__(self, feature_names):
        self.feature_names = feature_names
        self.scaler = StandardScaler()
        self.elastic_net1 = ElasticNet()
        self.elastic_net2 = ElasticNet()
        self.poly_features = None
        self.non_zero_coefs_ = None

    @staticmethod
    def rmse(y_true, y_pred):
        return np.sqrt(np.mean((y_pred - y_true) ** 2))

    @staticmethod
    def adjusted_r_squared(r_squared, n, p):
        return 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

    @staticmethod
    def mae(y_true, y_pred):
        return np.mean(np.abs(y_pred - y_true))

    def preprocess_data(self, df):
        X = df.dropna(subset=self.feature_names)
        y = X.pop('Age')
        X_scaled = self.scaler.fit_transform(X)
        return pd.DataFrame(X_scaled, columns=self.feature_names), y

    def fit_initial_elastic_net(self, X, y):
        self.elastic_net1.fit(X, y)
        self.non_zero_coefs_ = np.nonzero(self.elastic_net1.coef_)[0]

    def generate_polynomial_features(self, X):
        X_selected = X[:, self.non_zero_coefs_]
        self.poly_features = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
        X_poly = self.poly_features.fit_transform(X_selected)
        return X_poly

    def fit_second_elastic_net(self, X_poly, y):
        self.elastic_net2.fit(X_poly, y)
        # Update the non_zero_coefs_ to reflect the features selected by the second Elastic Net
        self.non_zero_coefs_ = np.nonzero(self.elastic_net2.coef_)[0]

    def train_model(self, X, y, alpha_range, l1_ratio_range, max_iter=100000):
        param_grid = {
            'alpha': alpha_range,
            'l1_ratio': l1_ratio_range,
            'max_iter': [max_iter]
        }
        rkfold = RepeatedKFold(n_splits=5, n_repeats=5)
        rmse_scorer = make_scorer(self.rmse, greater_is_better=False)
        gsearch = GridSearchCV(self.elastic_net2, param_grid, cv=rkfold, scoring=rmse_scorer, verbose=1, return_train_score=True)
        gsearch.fit(X, y)
        self.elastic_net2 = gsearch.best_estimator_

    def fit(self, df):
        X, y = self.preprocess_data(df)
        self.fit_initial_elastic_net(X, y)
        X_poly = self.generate_polynomial_features(X.values)
        self.fit_second_elastic_net(X_poly, y)
        # Here we could use a grid search as done in the train_model method
        # For simplicity, I'm calling the method directly without hyperparameter tuning
        alpha_range = np.logspace(-4, 2, 50)
        l1_ratio_range = np.arange(0.1, 1.0, 0.1)
        self.train_model(X_poly, y, alpha_range, l1_ratio_range)

    def predict(self, X):
        if self.non_zero_coefs_ is None or self.poly_features is None:
            raise AttributeError("The model has not been fitted yet. Call the fit method before predict.")
        X_scaled = self.scaler.transform(X)
        X_selected = X_scaled[:, self.non_zero_coefs_]
        X_poly = self.poly_features.transform(X_selected)
        return self.elastic_net2.predict(X_poly)

    def calculate_metrics(self, y_true, y_pred):
        r_squared = r2_score(y_true, y_pred)
        adjusted_r2 = self.adjusted_r_squared(r_squared, y_true.shape[0], len(self.non_zero_coefs_))
        rmse = self.rmse(y_true, y_pred)
        mae = self.mae(y_true, y_pred)
        return r_squared, adjusted_r2, rmse, mae

    def save_model(self, filename):
        # Implementation for saving the model to a file
        pass

    def load_model(self, filename):
        # Implementation for loading the model from a file
        pass

#######
predictor = iLGSScorePredictor(feature_names)
X, y = predictor.preprocess_data(df)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

predictor.fit(df)
y_pred = predictor.predict(X_test)
r_squared, rmse = predictor.calculate_metrics(y_test, y_pred)
predictor.plot_results(y_test, y_pred)

r_squared, rmse
