# predictSDEreversal
Attempting to predict the time until a reversal using machine learning techniques.

## Data

Data is previously computed in MATLAB. Data is an integration of an SDE of the form

<img src="https://latex.codecogs.com/gif.latex?dX_t=f(t)dt&space;&plus;&space;g(x)d\eta_t" title="dX_t=f(t)dt + g(x)d\eta_t" />

where the drift term is

<img src="https://latex.codecogs.com/gif.latex?f(x)&space;=&space;-\gamma*(|x|&space;-&space;x_0)sgn(x)(1-e^{-|x|/\epsilon})" title="f(x) = -\gamma*(|x| - x_0)sgn(x)(1-e^{-|x|/\epsilon})" />

the noise term is

<img src="https://latex.codecogs.com/gif.latex?g(x)&space;=&space;D" title="g(x) = D" />

and the noise process <img src="https://latex.codecogs.com/gif.latex?\eta" title="\eta" /> is an Ornstein-Uhlenbeck process with characteristic time <img src="https://latex.codecogs.com/gif.latex?\theta" title="\theta" />.

These values are fitted through comparison to the paleomagnetic model PADM2M ([Ziegler et al., 2011](https://doi.org/10.1111/j.1365-246X.2010.04905.x)). Drift and noise functions are displayed below.

<img src="/data/model.png" height="400"/>

Observed quantity is the amplitude of the axial dipole field. The quantity to be predicted from this is the "time until next polarity reversal".

<img src="/data/time.png" height="600"/>

Statistical features are calculated using a moving window over the ADM.

<img src="https://user-images.githubusercontent.com/38541020/90938056-032fa680-e3bd-11ea-902d-8b06acc60091.png" height="600"/>

## Predictions

A random forest method similar to the methodology of [Rouet‐Leduc et al., 2017](https://doi.org/10.1002/2017GL074677) is used. Hyperparameters have not been explored or tuned yet.

Preliminary results are not great. I believe this is because there is insufficient information in the SDE integration to accurately predict reversals. Physically richer models may provide the depth of information needed for predictions.

<img src="https://user-images.githubusercontent.com/38541020/90938118-2bb7a080-e3bd-11ea-9bee-dbea5e2f3657.png" height="600"/>


