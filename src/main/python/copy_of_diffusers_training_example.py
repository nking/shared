# -*- coding: utf-8 -*-
"""Copy of diffusers_training_example.ipynb
Automatically generated by Colaboratory.
Original file is located at
    https://colab.research.google.com/drive/1iXXHtz9lNPLXBLGyAVIXym_ab2aiCiHL

NLK:
* have edited to run a very small subset all the way through to edit for updated
arguments and data directories.  Have temporarily removed the cloud options, but
will enable them again.
* have added some math details to tie into background with noise filtering algorithms 
* have added a method for noise-filtering (non generative), that is, image restoration,
to explore whether the algorithm can remove noise from the image (like BM3D, etc.).
The authors also state at end of Sect 3.2 that when their DDPM model is configured to
predict the uncorrupted unknown image instead of predicting the noise, they found that
it worsened the sample quality.
The term "de-noise" used regarding their algorithm, seems to be used because the
DDPM loss objective resembles denoising score matching over multiple score scales indexed by t
as they state in Sect 3.2.

The DDPM paper reference is
Ho, J., Jain, A. and Abbeel, P., 2020. Denoising diffusion probabilistic models.
Advances in Neural Information Processing Systems, 33, pp.6840-6851.

> * This tutorial's choice of training images is a set of well curated images of butterflies.
One can also find a pretrained model from a set of cat images too "google/ddpm-cat-256".
The training uses many images that are carefully prepared to represent samples from the same distribution.
The DDPM paper authors (Ho, Jain, Abbeel 2020) however, use images that are busy, uncentered,
and do not have background removed.
This distribution is part of the encoding for the latent space.

B.T.W. The white background of the training images makes it easier to see effects of
applying the reverse process to the images (non-generative use case).

Some details on noise removal algorithms:

The form of "grouping and collaborative filtering" to learn noise and remove it can be seen in
the BM3D algorithm.
Others are wavelet transforms, Wiener, median, shrinkage-thresholding, linear and non-linear 
smoothing, statistical methods...

a simple form of deconvolution is (https://en.m.wikipedia.org/wiki/Deconvolution):
    f * g = h
       where f is source signal to recover,
       g is distortion function or filter,
       h is the observed signal
       * is the convolution symbol

    (f * g) + eps = h
       where eps is added noise

    if eps is small (or removed), we can recover f using fourier transforms
    followed by inverse fourier transform
       H = FFT(h); G = FFT(g)

    When modelling the noise to remove it, a common technique is to use a Weiner filter
    which assumes white noise, that is, a random signal with equal intensities at different frequencies
    (== constant power spectral density).

    DDPM uses a random signal and a random sampling of frequencies from a 
    bounded frequency domain in its estimation of noise to add to images.  
    This is part of the 'forward' process of DDPM 
    and is a cumulated product with re-parameterized parameters.

Ho, J., Jain, A. and Abbeel, P., 2020. Denoising diffusion probabilistic models.
Advances in Neural Information Processing Systems, 33, pp.6840-6851.
   The forward process of training takes an image and continuously adds gaussian noise to it.
(remember that any underlying noise distribution will resemble a gaussian distribution
as n increases, where n is the number of draws).
   The model learns the residual noise between the latest image and the latest image before
this round of adding noise.  The loss is calculated and the model parameters are updated
using the gradients of the loss. 
   The trained model can be used to "reverse" the process of adding noise to an image.
They call the process "de-noising" because they optimize a loss objective that resembles
de-noising score matching over multiple score scales indexed by t.

   A tangent to explore why the model doesn't work well as a noise-filter in image restoration:
   The noise model is composed of an aggregate of noise added to batches of images.
When the reverse process sees what it recognizes as noise in the input image,
it does not try to remove it, but instead, tries to restore the noise pixels to a state it
believes is before the noise was added.  The representation of the previous state for those
pixels may be very different than an individual image's pixels' nearest neighbors when the
individual starting input image is not pure noise.
   When the individual starting input image is pure noise, the restoration of pixels 
follows what the model has seen before essentially and generates an interesting image
from the implicit forward posterior distribution and UNet + noise scheduler model of the noise
in a decoding phase. see Sect 4.4 of paper.
   To explore the concept more, we could restrict the images to being binary, k = 2,
and having n_pixels = 8*4.  The number of possible images to create is k^n_pixels = 4294967296.
The training and denoising both have a cyclic complication.
It would be difficult to tell when to stop the restoration because each image - noise
would be another true image from the training dataset.  So one would want the training
set to be some subset of the 4294967296 which have properties like shapes being more than
1 adjacent pixel... could perform dilation on all 4294967296 images and keep only the unique subset,
and consider other properties like discard images of mostly all 0s or all 1s... Then the
universe of train and test images would be feasible to create... haven't thought this 
through thoroughly... would need to greatly reduce num_train_timesteps...

   The paper equations are here for convenience:

     reverse:
       eqn (1)
       p(x_{t-1} | x_{t}) = N( x_{t-1}; mu_{theta}(x_{t}, t), sigma_{theta}(x_{t}, t) )
     forward:
       eqn (2)
       q(x_{t} | x_{t-1}) = N( x_{t}; sqrt(1-beta_{t}) * x_{t-1}, beta_{t} * I)

     training performed using ELBO:
       eqn (3)
       L = E_{q}[-log(p(x_{T})) 
                 - sum over t >=1 to T (log(p_{theta}(x_{t-1}|x_{t}) / q(x_{t}|x_{t-1})))]

       alpha_{t} = 1 - beta_{t}

       alpha_mean_{t} = product over s = 1 to t (alpha_{s}

       eqn (4)
       q(x_{t}|x_{0}) = N(x_{t}; sqrt(alpha_mean{t}) * x_{0}, (1 - alpha_mean{t}) * I)

     rewritten using Kullbach-Leibler Divergence, and grouping of terms in 3 segments:
       eqn (5)
       L = E_{q}[ DKL(q(x_{T}|x_{0}) || p(x_{T})) <--L_{T}
                  + sum over t >=1 to T (
                    DKL(q(x_{t-1}|x_{t},x_{0}) || p(x_{t-1}|x_{t})) <--L_{t-1}
                    )
                  - log(p_{theta}(x_{0}|x_{1}))   <-- L_{0}
                ]
      
       eqn (6)
       q(x_{t-1}|x_{t},x_{0} = N( x_{t-1}; mu_est_{t}(x_{t}|x_{0}), beta_est_{t}*I )
       eqn (7)
          where 
            mu_est_{t}(x_{t}|x_{0}) = 
              (sqrt(1-alpha_mean_{t-1})*beta_{t} * x_{0}/(1-alpha_mean_{t})
              + sqrt(alpha_t)*(1-alpha_mean_{t}) * x_{t}/(1-alpha_mean_{t})
          where
            beta_est_{t} = (1-alpha_mean_{t-1}) * beta_{t} / (1-alpha_mean_{t})

       the DKLs can be calculated with closed form expressions in Rao-Blackwellized fashion.
              
       the authors fix beta_{t} of forward process, making L_{T} const during training (ignorable).
       The L_{t-1} term is simplified to eqn (12) in Sec 3.2
           L_{t-1} = E_x0_eps[((beta_{t}^2)/(2 * sigma_{t}^2 * alpha_{t}*(1-alpha_mean_{t-1})))
                       * ||(eps - eps_theta*(sqrt(alpha_mean_{t})*x_{0} + sqrt(1-alpha_mean_{t})*eps, t)||^2
                     ]

           where eps_theta is a function approximator to predict eps from x_{t}, where x_{t} is a noisey image.

       Sec 3.2 of paper:
       the authors set (sigma_{t})^2 to beta_{t} for generative use, x_{0}~N(0,I)
       ** but for other uses, could set (sigma_{t})^2 = beta_est_{t}

       in Sect 3.2 last paragraph, the authors state that they train the reverse process
       mean function approximator mu_theta to predict the noise, eps. 
       (they add that one can predict x_{0}, but found it led to worse sample quality).

       UNet model, specifically, is used for the mean function approximator mu_theta (which is
       configured to predict the noisy image).
       DDPMScheduler holds the alphas indexed by t.
       The training loss (empirical risk) is calculated as square of difference between noisey
       image created in reverse process and the noise predicted by UNet model.
       The gradient of the loss is performed by the optimizer, which updates the model parameters
       (where parameters are the weights and biases in the Model network layers).

#######
To run the code all the way though to check arguments for libraries you have installed,
edit run_small_inspect to be True below.
#######
# 🤗 Training with Diffusers

In recent months, it has become clear that diffusion models have taken the throne as the state-of-the-art generative models. Here, we will use Hugging Face's brand new [Diffusers](https://github.com/huggingface/diffusers) library to train a simple diffusion model.

## Installing the dependencies

This notebook leverages the [🤗 Datasets](https://huggingface.co/docs/datasets/index) library to load and preprocess image datasets and the [🤗 Accelerate](https://huggingface.co/docs/accelerate/index) library to simplify training on any number of GPUs, with features like automatic gradient accumulation and tensorboard logging. Let's install them here:
"""

# Commented out IPython magic to ensure Python compatibility.
# %%capture
# !pip install git+https://github.com/huggingface/diffusers.git#egg=diffusers[training]
# !pip install accelerate
# !pip install datasets

"""If you're opening this notebook locally, make sure your environment has an install from the last version of those libraries.

To be able to share your model with the community, there are a few more steps to follow.|

First you have to store your authentication token from the Hugging Face website (sign up [here](https://huggingface.co/join) if you haven't already!) then execute the following cell and input your **write** token:
"""

#from huggingface_hub import notebook_login

#notebook_login()

"""
Then you need to install Git-LFS to upload your model checkpoints:"""

# Commented out IPython magic to ensure Python compatibility.
# %%capture
# !sudo apt -qq install git-lfs
# !git config --global credential.helper store

"""## Config

For convenience, we define a configuration grouping all the training hyperparameters. This would be similar to the arguments used for a [training script](https://github.com/huggingface/diffusers/tree/main/examples).
Here we choose reasonable defaults for hyperparameters like 
`num_epochs`, `learning_rate`, `lr_warmup_steps`, 
but feel free to adjust them if you train on your own dataset. 
For example, `num_epochs` can be increased to 100 for better visual quality.
"""

# to run a small subset of the data, small number of epochs and iterations
# in order to follow the process edit data directories, etc set this to True:
run_small_inspect = True

if run_small_inspect:
    print(f'WARNING: run_small_inspect=True, (uses small subset of data, small # train epochs small # sample iterations')

display_train_sample = False


import os

hf_name = "huggan/smithsonian_butterflies_subset"

try:
    data_dir = os.environ["ML_DATASETS_HOME"]
    out_data_dir = os.environ["ML_OUTPUT_HOME"]
except:
    import tempfile
    from tempfile import TemporaryDirectory
    # colab environment
    # other args: url, folder_in_archive
    data_dir = tempfile.mkdtemp()
    print(f'temp_dataset_dir={data_dir}')
    out_data_dir = os.path.join(data_dir, "output")

data_dir = os.path.join(data_dir, hf_name)
out_data_dir = os.path.join(out_data_dir, hf_name)
if run_small_inspect:
    out_data_dir = out_data_dir + "_small"

print(f'data_dir is {data_dir}')
print(f'out_data_dir is {out_data_dir}')

from dataclasses import dataclass

@dataclass
class TrainingConfig:
    image_size = 128  # the generated image resolution
    # have reduced the batch sizes to reduce use of RAM
    train_batch_size = 4#16
    eval_batch_size = 4#16  # how many images to sample during evaluation
    num_epochs = 50
    gradient_accumulation_steps = 1
    learning_rate = 1e-4
    lr_warmup_steps = 500
    save_image_epochs = 10
    save_model_epochs = 30
    mixed_precision = 'fp16'  # `no` for float32, `fp16` for automatic mixed precision
    dataset = hf_name
    output_dir = out_data_dir  # the model name locally and on the HF Hub

    push_to_hub = False  # whether to upload the saved model to the HF Hub
    hub_private_repo = False
    overwrite_output_dir = True  # overwrite the old model when re-running the notebook
    seed = 0
    num_train_timesteps = 1000

config = TrainingConfig()

if run_small_inspect:
    config.num_epochs = 5
    config.train_batch_size = 4
    config.train_batch_size = 4
    config.eval_batch_size = 4
    config.num_train_timesteps = 10

"""## Loading the dataset

We will use the [🤗 Datasets](https://github.com/huggingface/datasets) library to download our image dataset.

In this case, the [Butterflies dataset](https://huggingface.co/datasets/huggan/smithsonian_butterflies_subset) is hosted remotely, but you can load a local [ImageFolder](https://huggingface.co/docs/datasets/v2.0.0/en/image_process#imagefolder) as shown in the commets below.
"""
from datasets import load_dataset

if (not os.path.exists(data_dir)):
    # butterflies
    dataset = load_dataset(path=config.dataset, split="train", keep_in_memory=False)
    dataset.save_to_disk(data_dir)
else:
    from datasets import Dataset
    #from torch.utils.data import Dataset
    # keep_in_memory = False?
    dataset = Dataset.load_from_disk(data_dir, keep_in_memory=False)

if run_small_inspect:
    # load only 1 batch to check pipeline
    import torch
    import numpy as np
    ind = np.array([i for i in range(0, config.train_batch_size)], dtype=np.int32)
    indices = torch.from_numpy(ind)
    ind = ind + config.train_batch_size
    indices2 = torch.from_numpy(ind)
    dataset_train = torch.utils.data.Subset(dataset, indices)
    #dataset_train_ = torch.utils.data.RandomSampler(dataset, num_samples=config.train_batch_size)
    #dataset_test = torch.utils.data.Subset(dataset, indices2)
    dataset = dataset_train

# Feel free to try other datasets from https://hf.co/huggan/ too!
# Here's is a dataset of flower photos:
# config.dataset = "huggan/flowers-102-categories"
# dataset = load_dataset(config.dataset, split="train")

'''
The dataset contains several extra `features` (columns), but the one that we're interested in is `image`:
Dataset({
    features: ['image_url', 'image_alt', 'id', 'name', 'scientific_name', 'gender', 'taxonomy', 'region', 'locality', 
    'date', 'usnm_no', 'guid', 'edan_url', 'source', 'stage', 'image', 'image_hash', 'sim_score'],
    num_rows: 1000
})
'''

"""Since the [`Image`](https://huggingface.co/docs/datasets/image_process#image-datasets) feature loads the images with PIL, we can easily look at a few examples:"""

import matplotlib.pyplot as plt

if display_train_sample:
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    for i, image in enumerate(dataset[:4]["image"]):
        axs[i].imshow(image)
        axs[i].set_axis_off()
    fig.show()
else:
    print('skipping plot figures this time')

"""The images in the dataset are all different, so we need to preprocess them first:
* `Resize` makes the images conform to a square resolution of `config.image_size`
* `RandomHorizontalFlip` augments the dataset by randomly mirroring the images.
* `Normalize` is important to rescale the pixel values into a `[-1, 1]` range (which our model will expect).
"""
#!pip3 install torchvision
from torchvision import transforms

preprocess = transforms.Compose(
    [
        transforms.Resize((config.image_size, config.image_size)),
        transforms.RandomHorizontalFlip(),
        transforms.ToTensor(),
        transforms.Normalize([0.5], [0.5]),
    ]
)

"""🤗 Datasets offer a handy `set_transform()` method to apply the image transformations on the fly during training:"""

def transform(examples):
    images = [preprocess(image.convert("RGB")) for image in examples["image"]]
    return {"images": images}

try:
    dataset.set_transform(transform)
except:
    #for Subset
    dataset.dataset.set_transform(transform)

"""Let's see what they look like now"""

if display_train_sample:
    fig, axs = plt.subplots(1, 4, figsize=(16, 4))
    for i, image in enumerate(dataset[:4]["images"]):
        axs[i].imshow(image.permute(1, 2, 0).numpy() / 2 + 0.5)
        axs[i].set_axis_off()
    fig.show()
else:
    print('skipping plot figures this time')

"""Now that all our images have the same size and are converted to tensors, we can create the dataloader we will use for training."""

import torch

train_dataloader = torch.utils.data.DataLoader(dataset, batch_size=config.train_batch_size, shuffle=True)

"""## Defining the diffusion model

Here we set up our diffusion model. Diffusion models are neural networks that are trained to predict slightly less noisy images from a noisy input. At inference, they can be used to iteratively transform a random noise to generate an image:

<p align="center">
    <img src="https://user-images.githubusercontent.com/10695622/174349667-04e9e485-793b-429a-affe-096e8199ad5b.png" width="800"/>
    <br>
    <em> Figure from DDPM paper (https://arxiv.org/abs/2006.11239). </em>
<p>

Don't worry too much about the math if you're not familiar with it, the import part to remember is that our model corresponds to the arrow $p_{\theta}(x_{t-1}|x_{t})$ (which is a fancy way of saying: predict a slightly less noisy image).

The interesting part is that it's really easy to add some noise to an image, so the training can happen in a semi-supervised fashion as follows:
1. Take an image from the training set.
2. Apply to it some random noise $t$ times (this will give the $x_{t-1}$ and the $x_{t}$ in the figure above).
3. Give this noisy image to the model along with the value of $t$.
4. Compute a loss from the output of the model and the noised image $x_{t-1}$.

Then we can apply gradient descent and repeat this process multiple times.

Most diffusion models use architectures that are some variant of a [U-net](https://arxiv.org/abs/1505.04597) and that's what we'll use here.

![](https://huggingface.co/datasets/huggingface/documentation-images/resolve/main/unet-model.png)

In a nutshell:
- the model has the input image go through several blocks of ResNet layers which halves the image size by 2
- then through the same number of blocks that upsample it again.
- there are skip connections linking the features on the downample path to the corresponsding layers in the upsample path.

A key feature of this model is that it predicts images of the same size as the input, which is exactly what we need here.

Diffusers provides us a handy `UNet2DModel` class which creates the desired architecture in PyTorch.

Let's create a U-net for our desired image size.
Note that `down_block_types` correspond to the downsampling blocks (green on the diagram above), and `up_block_types` are the upsampling blocks (red on the diagram):
"""

from diffusers import UNet2DModel

# documentation https://huggingface.co/docs/diffusers/api/models/unet2d
model = UNet2DModel(
    sample_size=config.image_size,  # the target image resolution
    in_channels=3,  # the number of input channels, 3 for RGB images
    out_channels=3,  # the number of output channels
    layers_per_block=2,  # how many ResNet layers to use per UNet block
    block_out_channels=(128, 128, 256, 256, 512, 512),  # the number of output channes for each UNet block
    down_block_types=(
        "DownBlock2D",  # a regular ResNet downsampling block
        "DownBlock2D",
        "DownBlock2D",
        "DownBlock2D",
        "AttnDownBlock2D",  # a ResNet downsampling block with spatial self-attention
        "DownBlock2D",
    ),
    up_block_types=(
        "UpBlock2D",  # a regular ResNet upsampling block
        "AttnUpBlock2D",  # a ResNet upsampling block with spatial self-attention
        "UpBlock2D",
        "UpBlock2D",
        "UpBlock2D",
        "UpBlock2D"
      ),
)

"""Let's get a sample image from our dataset and pass it into our model. We just need to add a batch dimension:"""

try:
    sample_image = dataset[0]['images'].unsqueeze(0)
except:
    #Subset
    sample_image = dataset.dataset[0]['images'].unsqueeze(0)
print('Input shape:', sample_image.shape)

"""And let's check the output is a tensor of the same exact shape:"""

print('Output shape:', model(sample_image, timestep=0)["sample"].shape, flush=True)

"""Great!

Note that our model takes in the (noisy) image and also the current time-step (as we saw before in the training overview). That time-step information is converted for the model using a sinusoidal positional embedding, similar to what Transformer models often do.

Now that we have our model, we just need an object to *add noise to an image*. This is done by the **schedulers** in the `diffusers` library.

## Defining the noise scheduler

Depending on the diffusion algorithm you want to use, the way images are noised is slightly different. That's why 🤗 Diffusers contains different scheduler classes which each define the algorithm-specific diffusion steps. Here we are going to use the `DDPMScheduler` which corresponds to the training denoising and training algorithm proposed in [Denoising Diffusion Probabilistic Models](https://arxiv.org/abs/2006.11239).
"""

from diffusers import DDPMScheduler

noise_scheduler = DDPMScheduler()
noise_scheduler.config.num_train_timesteps = config.num_train_timesteps
'''
other args for DDPMScheduler:
prediction_type (`str`, defaults to `epsilon`, *optional*):
            Prediction type of the scheduler function; can be `epsilon` (predicts the noise of the diffusion process),
            `sample` (directly predicts the noisy sample`) or `v_prediction` (see section 2.4 of [Imagen
            Video](https://imagen.research.google/video/paper.pdf) paper).
            
'''

"""Let's see how this noise scheduler works: it takes a batch of images from the trainng set (here we will reuse the 
batch of one image `sample_image` form before), a batch of random noise of the same shape 
and the timesteps for each image (which correspond to the number of times we want to apply noise to each image):"""

import torch
from PIL import Image

noise = torch.randn(sample_image.shape)
timesteps = torch.LongTensor([50])
noisy_image = noise_scheduler.add_noise(sample_image, noise, timesteps)

Image.fromarray(((noisy_image.permute(0, 2, 3, 1) + 1.0) * 127.5).type(torch.uint8).numpy()[0])

"""In the DDPM algorithm, the training objective of the model is then to be able to predict the noise we used in `noise_scheduler.add_noise`, so the loss at this step would be:"""

import torch.nn.functional as F

noise_pred = model(noisy_image, timesteps)["sample"]
loss = F.mse_loss(noise_pred, noise)

"""## Setting up training

We have all we need to be able to train our model! Let's use a standard AdamW optimizer:
"""

optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate)

""" And a cosine learning rate schedule:"""

from diffusers.optimization import get_cosine_schedule_with_warmup

if run_small_inspect:
    # I change train loop to change lr step only once per epoch
    lr_scheduler = get_cosine_schedule_with_warmup(
        optimizer=optimizer,
        num_warmup_steps=config.lr_warmup_steps,
        num_training_steps=(config.num_epochs),
    )
else:
    #TODO: edit this when re-enable use of Accelerator
    lr_scheduler = get_cosine_schedule_with_warmup(
        optimizer=optimizer,
        num_warmup_steps=config.lr_warmup_steps,
        #num_training_steps=(len(train_dataloader) * config.num_epochs),
        num_training_steps=(config.num_epochs),
    )

"""To evaluate our model, we use the `DDPMPipeline` which is an easy way to perform end-to-end inference (see this notebook [TODO link] for more detail). We will use this pipeline to generate a batch of sample images and save it as a grid to the disk."""

from diffusers import DDPMPipeline

import math

def make_grid(images, rows, cols):
    w, h = images[0].size
    grid = Image.new('RGB', size=(cols*w, rows*h))
    for i, image in enumerate(images):
        grid.paste(image, box=(i%cols*w, i//cols*h))
    return grid

def evaluate(config, epoch, pipeline):
    # Sample some images from random noise (this is the backward diffusion process).
    # The default pipeline output type is `List[PIL.Image]`
    print(f'evaluate {epoch}')
    # default num_inference_steps=1000

    result = pipeline(
        batch_size = config.eval_batch_size,
        generator=torch.manual_seed(config.seed),
        #TODO: should use pipeline scheduler's config.num_train_timesteps 
        num_inference_steps=noise_scheduler.config.num_train_timesteps,
    )

    # result type is ImagePipelineOutput
    #print(f'result={result}')
    images = result["images"]

    # Make a grid out of the images
    image_grid = make_grid(images, rows=4, cols=4)

    # Save the images
    print('save sampled images to disk')
    test_dir = os.path.join(config.output_dir, "samples")
    os.makedirs(test_dir, exist_ok=True)
    image_grid.save(f"{test_dir}/{epoch:04d}.png")

"""With this in end, we can group all together and write our training function. 
This just wraps the training step we saw in the previous section in a loop, 
using Accelerate for easy TensorBoard logging, gradient accumulation, 
mixed precision training and multi-GPUs or TPU training."""

# https://pytorch.org/tutorials/intermediate/tensorboard_tutorial.html

from accelerate import Accelerator
import diffusers
#from diffusers.utils.hub_utils import init_git_repo, push_to_hub

from tqdm.auto import tqdm
import os
'''
def train_loop(config, model, noise_scheduler, optimizer, train_dataloader, lr_scheduler):
    # Initialize accelerator and tensorboard logging
    accelerator = Accelerator(
        mixed_precision=config.mixed_precision,
        gradient_accumulation_steps=config.gradient_accumulation_steps,
        log_with="tensorboard",
        #logging_dir=os.path.join(config.output_dir, "logs")
        project_dir=os.path.join(config.output_dir, "logs")
    )
    if accelerator.is_main_process:
        #if config.push_to_hub:
        #    repo = init_git_repo(config, at_init=True)
        accelerator.init_trackers("train_example")

    # Prepare everything
    # There is no specific order to remember, you just need to unpack the
    # objects in the same order you gave them to the prepare method.
    model, optimizer, train_dataloader, lr_scheduler = accelerator.prepare(
        model, optimizer, train_dataloader, lr_scheduler
    )

    global_step = 0

    # Now you train the model
    for epoch in range(config.num_epochs):
        progress_bar = tqdm(total=len(train_dataloader), disable=not accelerator.is_local_main_process)
        progress_bar.set_description(f"Epoch {epoch}")

        for step, batch in enumerate(train_dataloader):
            clean_images = batch['images']
            # Sample noise to add to the images
            noise = torch.randn(clean_images.shape).to(clean_images.device)
            bs = clean_images.shape[0]

            # Sample a random timestep for each image
            timesteps = torch.randint(0, noise_scheduler.num_train_timesteps, (bs,), device=clean_images.device).long()

            # Add noise to the clean images according to the noise magnitude at each timestep
            # (this is the forward diffusion process)
            noisy_images = noise_scheduler.add_noise(clean_images, noise, timesteps)

            with accelerator.accumulate(model):
                # Predict the noise residual
                noise_pred = model(noisy_images, timesteps)["sample"]
                loss = F.mse_loss(noise_pred, noise)
                accelerator.backward(loss)

                accelerator.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
                lr_scheduler.step()
                optimizer.zero_grad()

            progress_bar.update(1)
            logs = {"loss": loss.detach().item(), "lr": lr_scheduler.get_last_lr()[0], "step": global_step}
            progress_bar.set_postfix(**logs)
            accelerator.log(logs, step=global_step)
            global_step += 1

        # After each epoch you optionally sample some demo images with evaluate() and save the model
        if accelerator.is_main_process:
            pipeline = DDPMPipeline(unet=accelerator.unwrap_model(model), scheduler=noise_scheduler)

            if (epoch + 1) % config.save_image_epochs == 0 or epoch == config.num_epochs - 1:
                evaluate(config, epoch, pipeline)

            if (epoch + 1) % config.save_model_epochs == 0 or epoch == config.num_epochs - 1:
                #if config.push_to_hub:
                #    push_to_hub(config, pipeline, repo, commit_message=f"Epoch {epoch}", blocking=True)
                #else:
                #    pipeline.save_pretrained(config.output_dir)
                pipeline.save_pretrained(config.output_dir)
'''

def train_loop_no_accelerator(config, model, noise_scheduler, optimizer, train_dataloader, lr_scheduler):
    print('start train loop')
    #config.mixed_precision,
    #config.gradient_accumulation_steps,

    size = len(train_dataloader.dataset)
    model.train()  # sets the module to training mode

    # Now you train the model
    for epoch in range(config.num_epochs):

        for step, batch in enumerate(train_dataloader):
            clean_images = batch['images'] #(16,3,128,128)
            # Sample noise to add to the images
            noise = torch.randn(clean_images.shape).to(clean_images.device)
            bs = clean_images.shape[0]

            # Sample a random timestep for each image
            timesteps = torch.randint(0, noise_scheduler.config.num_train_timesteps, (bs,), device=clean_images.device).long()

            # Add noise to the clean images according to the noise magnitude at each timestep
            # (this is the forward diffusion process)
            # alphas are kept in noise_scheduler.
            noisy_images = noise_scheduler.add_noise(clean_images, noise, timesteps)

            # Predict the noise residual
            noise_pred = model(noisy_images, timesteps)["sample"]
            loss = F.mse_loss(noise_pred, noise)

            # TODO: consider plotting (clean_images - noise_pred) as training progresses.

            # Backpropagation
            if run_small_inspect:
                print('backpropagation', flush=True)
            loss.backward()
            if run_small_inspect:
                print('update params', flush=True)
            # update the parameters (weights and biases)  of UNet2DModel indirectly through optimizer
            optimizer.step()
            optimizer.zero_grad()

            if run_small_inspect or step % 100 == 0:
                loss = loss.detach().item()
                print(f"loss: {loss:>7f}  [{step:>5d}/{size:>5d}, {epoch}]")

        # after an epoch of train and eval, let the regularization adapt:
        lr_scheduler.step()

        # if this is the last epoch, run the computationally expensive evaluation (sampling is expensive).
        if (epoch + 1) % config.save_image_epochs == 0 or epoch == config.num_epochs - 1:
            print('construct pipeline')
            # pipeline = DDPMPipeline(unet=accelerator.unwrap_model(model), scheduler=noise_scheduler)
            pipeline = DDPMPipeline(unet=model, scheduler=noise_scheduler)

            # moved evaluate to outside of epochs loop

            # if config.push_to_hub:
            #    push_to_hub(config, pipeline, repo, commit_message=f"Epoch {epoch}", blocking=True)
            # else:
            #    pipeline.save_pretrained(config.output_dir)
            print('save pipeline configuration')
            pipeline.save_pretrained(config.output_dir)

"""## Let's train!

Let's launch the training (including multi-GPU training) from the notebook using Accelerate's `notebook_launcher` function:
"""

from accelerate import notebook_launcher
args = (config, model, noise_scheduler, optimizer, train_dataloader, lr_scheduler)

#notebook_launcher(train_loop, args, num_processes=1)
train_loop_no_accelerator(config, model, noise_scheduler, optimizer, train_dataloader, lr_scheduler)
print('done training')

import glob

def generate_sample() :
    '''
    apply the trained U-Net model, UNet2DModel, and alpha coefficients to images of 
    pure noise.  the result is the creation of signal in the image. the final image
    belongs to same probability distribution as the original images.
    '''
    print('evaluate generative pipeline')
    pipeline = DDPMPipeline(unet=model, scheduler=noise_scheduler)
    evaluate(config, config.num_epochs, pipeline)
    sample_images = sorted(glob.glob(f"{config.output_dir}/samples/*.png"))
    Image.open(sample_images[-1])
    print('done with display sample generated images')

def filter_noise():
    import math
    '''
    apply the trained U-Net model, UNet2DModel, and alpha coefficients to the original images
    to restore an image corrupted by noise.  The input images used are the training images.
    one could substitute other images here.

    TODO: the images possibly need to be pre-processed to coordinate-wise unit variance.
          they've been scaled by a normalize operation only.
    TODO: the sigma_t term possibly needs to be set to beta_t * (1- (alpha_bar_{t-1}))/(1-(alpha_bar_{t}))
          see Ho et al. 2020 eqn 8.
    '''
    dirname = "de-noise"
    img_dir = os.path.join(config.output_dir, dirname)
    os.makedirs(img_dir, exist_ok=True)
    print(f'apply de-noise filtering to original images')

    timesteps = torch.from_numpy(np.arange(0, noise_scheduler.config.num_train_timesteps)[::-1].copy())

    for step, batch in enumerate(train_dataloader):
        image_batch = batch['images']  # (16,3,128,128)

        ## diffusers.pipelines.pipeline_utils.DiffusionPipeline.progress_bar
        for t in tqdm(timesteps):
            # 1. predict noise model_output
            model_output = model(image_batch, t).sample

            # 2. compute previous image: x_t -> x_t-1 by reversing the stochastic differential equation.
            # This function propagates the diffusion process from the learned model outputs
            # (most often the predicted noise).
            image_batch = noise_scheduler.step(model_output, t, image_batch,
                    generator=torch.manual_seed(config.seed),).prev_sample

            if t % 100 == 0:
                img = (image_batch / 2 + 0.5).clamp(0, 1)
                img = img.cpu().permute(0, 2, 3, 1).detach().numpy()
                img = diffusers.pipelines.pipeline_utils.numpy_to_pil(img)
                image_grid = make_grid(img, rows=4, cols=4)
                image_grid.save(f"{img_dir}/{step:04d}_t{t:04d}.png")
                print(f"wrote to {img_dir}/{step:04d}_t{t:04d}.png")

        image_batch = (image_batch / 2 + 0.5).clamp(0, 1)
        image_batch = image_batch.cpu().permute(0, 2, 3, 1).detach().numpy()
        image_batch = diffusers.pipelines.pipeline_utils.numpy_to_pil(image_batch)

        # Make a grid out of the images
        image_grid = make_grid(image_batch, rows=4, cols=4)

        # Save the images
        print('save noise-filtered images {step} to disk')
        image_grid.save(f"{img_dir}/{step:04d}.png")
        print(f"wrote to {img_dir}/{step:04d}.png")

    filtered_images = sorted(glob.glob(f"{img_dir}/*.png"))
    #Image.open(filtered_images[-1])

    nt = len(filtered_images)
    ncols = 1
    nrows = int(math.ceil(nt/ncols))
    j = 1
    for idx in range(0, nt, 1):
        img = np.asarray(Image.open(filtered_images[idx]))
        ax = plt.subplot(nrows, ncols, j)
        plt.imshow(img)
        j += 1


    plt.show()
    print('done with display noise-filtered images')

# these are 2 seperate use-cases for the trained model
filter_noise()
generate_sample()

if run_small_inspect:
    print(f'WARNING REMINDER: run_small_inspect=True, (uses small subset of data, small # train epochs small # sample iterations')

"""Not bad! There's room for improvement of course, so feel free to play with the hyperparameters, model definition and image augmentations 🤗

If you've chosen to upload the model to the Hugging Face Hub, its repository should now look like so:
https://huggingface.co/anton-l/ddpm-butterflies-128

If you want to dive deeper into the code, we also have more advanced training scripts with features like Exponential Moving Average of model weights here:

https://github.com/huggingface/diffusers/tree/main/examples
"""

