# AuDioMio: Audio BPM Estimation Tool

This repository contains a simple Rust tool for estimating the Beats Per
Minute (BPM) of an audio file. The tool processes the audio file, calculates its
spectrogram, detects peaks, and utilizes histogram analysis to provide a BPM 
estimation.

I built it to learn a little bit about digital signal processing.

## Introduction

Estimating BPM (Beats Per Minute) in music is a fundamental task in audio
processing and analysis. This tool uses Rust to provide a basic yet effective
implementation for estimating BPM from audio files.

## Algorithms Used

### Spectrogram Calculation

#### What is a Spectrogram?

A spectrogram is a visual representation of the frequencies present in a signal
over time. In the context of audio, it helps us see how different frequencies
contribute to the overall sound. Think of it as a musical landscape where the
x-axis represents time, the y-axis represents frequency, and the color or
intensity represents the amplitude (loudness) of each frequency at a given
moment.

#### How is it Calculated?

To create a spectrogram, we use a mathematical tool called the Fast Fourier
Transform (FFT). The FFT takes a snapshot of the audio at small intervals,
converts it from the time domain to the frequency domain, and helps identify
which frequencies are present in that snapshot.

Additionally, we use a Hamming window to smooth out the edges of these
snapshots. The windowing process helps prevent a phenomenon known as spectral
leakage, which can introduce artifacts into the spectrogram.

### Peak Detection

#### What are Peaks?

Peaks in the context of audio represent moments where certain frequencies are
particularly strong or pronounced. In our case, we're interested in detecting
peaks in the spectrogram that correspond to beats or rhythmic elements in the
music.

#### How are Peaks Detected?

We apply a simple thresholding approach to identify peaks. If the amplitude of a
particular frequency at a given time exceeds a certain threshold, we consider it
a peak. These peaks help us locate significant moments in the audio where beats
are likely to occur.

### BPM Estimation with Histogram Analysis

#### What is a Histogram?

A histogram is like a bar chart that shows the distribution of values within a
certain range. In our case, we're interested in the distribution of time
intervals between detected peaks.

#### How is BPM Estimated?

We calculate the time intervals between consecutive peaks and create a histogram
to see which intervals occur most frequently. The idea is that the most common
interval likely corresponds to the underlying tempo of the music.

By analyzing this histogram, we can identify the dominant tempo and calculate
the Beats Per Minute (BPM) based on the average time interval between peaks.
This approach provides a robust estimation even if the tempo of the music
varies.

These algorithms work together to analyze the audio and provide an estimate of
the BPM, helping us understand the rhythm and tempo of the music.

## Usage

1. **Install Rust:**

   Ensure that you have Rust installed on your system. You can install it by
   following the instructions on
   the [official Rust website](https://www.rust-lang.org/tools/install).

2. **Clone the Repository:**

   ```bash
   git clone https://github.com/simpajj/audio-bpm-estimation.git
   cd audio-bpm-estimation
   ```

3. **Build and Run:**

   ```bash
   cargo build --release
   cargo run --release -- <path_to_your_audio_file.wav>
   ```

   > Replace <path_to_your_audio_file.wav> with the actual path to your audio 
   > file in WAV format.

The tool will process the audio file and output the estimated BPM.
