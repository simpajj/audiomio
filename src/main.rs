use std::env;
use hound::Error;
use rustfft::FftPlanner;
use rustfft::num_complex::Complex;

// The size of the analysis window.
const WINDOW_SIZE: usize = 1024;
// The step size between consecutive windows.
const HOP_SIZE: usize = WINDOW_SIZE / 2;
// The threshold for peak detection.
const THRESHOLD: f32 = 0.1;

struct AudioTrack {
  pub samples: Vec<f32>,
  pub sample_rate: u32,
}

fn main() -> Result<(), Error> {
  let file_path = env::args().nth(1).expect("No file path given");
  let track = read_audio_file(&file_path)?;

  let spectrogram = calculate_spectrogram(&track.samples);
  let peaks = detect_peaks(&spectrogram);
  let bpm = calculate_bpm(&peaks, track.sample_rate);

  println!("Estimated BPM: {:.2}", bpm);

  Ok(())
}

/// Reads an audio file and returns a vector of f32 samples.
///
/// Reading the audio file is the initial step in the analysis process.
/// It provides the raw data needed for further processing.
///
/// # Arguments
///
/// * `file_path` - The path to the audio file in WAV format.
///
/// # Returns
///
/// A Result containing an AudioTrack with a vector of f32 samples
/// and the sample rate if successful, or an error if reading the
/// file fails.
fn read_audio_file(file_path: &str) -> Result<AudioTrack, Error> {
  let mut reader = hound::WavReader::open(file_path)?;
  let samples: Vec<f32> = reader.samples::<i16>().map(|s| s.unwrap() as f32).collect();
  Ok(AudioTrack {
    samples,
    sample_rate: reader.spec().sample_rate,
  })
}

/// Calculates the spectrogram of an audio signal using FFT
/// and a Hamming window.
///
/// The spectrogram provides a time-frequency representation
/// of the audio signal, essential for identifying rhythmic
/// patterns and beats.
///
/// # Arguments
///
/// * `samples` - A slice of audio samples.
///
/// # Returns
///
/// A 2D vector representing the magnitude spectrogram.
fn calculate_spectrogram(samples: &[f32]) -> Vec<Vec<f32>> {
  let fft_size = WINDOW_SIZE.next_power_of_two();
  let fft = FftPlanner::new().plan_fft_forward(fft_size);

  let mut spectrogram = Vec::new();

  for window_start in (0..samples.len()).step_by(HOP_SIZE) {
    let window_end = window_start + WINDOW_SIZE;
    if window_end > samples.len() {
      break;
    }

    let mut windowed_samples: Vec<f32> = samples[window_start..window_end].to_vec();
    normalize_audio_samples(&mut windowed_samples);

    let mut complex_samples: Vec<_> = windowed_samples.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fft.process(&mut complex_samples);

    let magnitude_spectrum: Vec<_> = complex_samples.iter().map(|&c| c.norm()).collect();
    spectrogram.push(magnitude_spectrum);
  }

  spectrogram
}

/// Applies the Hamming window to a vector of audio samples.
///
/// Applying a window function, such as Hamming, helps mitigate
/// spectral leakage during the Fourier Transform, improving the
/// accuracy of frequency analysis.
///
/// # Arguments
///
/// * `samples` - A mutable reference to a vector of audio samples.
fn normalize_audio_samples(samples: &mut Vec<f32>) {
  let window: Vec<f32> = (0..samples.len())
    .map(|i| 0.54 - 0.46 * f32::cos(2.0 * std::f32::consts::PI * i as f32 / (samples.len() - 1) as f32))
    .collect();

  for i in 0..samples.len() {
    samples[i] *= window[i];
  }
}

/// Detects peaks in a spectrogram using a simple thresholding approach.
///
/// Detecting peaks helps identify prominent features in the spectrogram,
///
/// # Arguments
///
/// * `spectrogram` - The magnitude spectrogram of an audio signal.
/// * `threshold` - The threshold for peak detection.
///
/// # Returns
///
/// A vector of tuples representing the coordinates of detected peaks (row, column) in the spectrogram.
fn detect_peaks(spectrogram: &[Vec<f32>]) -> Vec<(usize, usize)> {
  let mut peaks = Vec::new();

  for (i, row) in spectrogram.iter().enumerate() {
    for (j, &magnitude) in row.iter().enumerate() {
      if magnitude > THRESHOLD && is_local_maximum(spectrogram, i, j) {
        peaks.push((i, j));
      }
    }
  }

  peaks
}

/// Checks if a point in a spectrogram is a local maximum.
///
/// Identifying local maxima helps filter out insignificant
/// peaks and focus on prominent features in the spectrogram,
/// improving the accuracy of beat detection.
///
/// # Arguments
///
/// * `spectrogram` - The magnitude spectrogram of an audio signal.
/// * `i` - The row index of the point.
/// * `j` - The column index of the point.
///
/// # Returns
///
/// A boolean indicating whether the point is a local maximum.
fn is_local_maximum(spectrogram: &[Vec<f32>], i: usize, j: usize) -> bool {
  let magnitude = spectrogram[i][j];

  // Check neighbors for lower magnitude
  let neighbors = [
    (i.wrapping_sub(1), j),
    (i.wrapping_add(1), j),
    (i, j.wrapping_sub(1)),
    (i, j.wrapping_add(1)),
  ];

  // Check if the magnitude of the current point is greater than the magnitudes of all its neighbors
  neighbors.iter().all(|&(ni, nj)| {
    ni < spectrogram.len() && nj < spectrogram[i].len() && spectrogram[ni][nj] < magnitude
  })
}

/// Creates a histogram of time intervals.
fn create_histogram(intervals: &[f32], bins: usize) -> Vec<usize> {
  let min_interval = *intervals.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&0.0);
  let max_interval = *intervals.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap_or(&0.0);

  let bin_width = (max_interval - min_interval) / bins as f32;
  let mut histogram = vec![0; bins];

  for &interval in intervals {
    let bin_index = ((interval - min_interval) / bin_width).floor() as usize;

    // Ensure the bin_index is within bounds
    let bin_index = bin_index.min(bins - 1);

    histogram[bin_index] += 1;
  }

  histogram
}

/// Calculates the Beats Per Minute (BPM) from detected peaks.
///
/// Estimating BPM from detected peaks provides a quantitative
/// measure of the tempo of the audio. The function uses histogram
/// analysis to identify the most common time interval between peaks,
/// providing a more robust and accurate estimation in cases where tempo
/// variations exist within the audio.
///
/// # Arguments
///
/// * `peaks` - A vector of tuples representing the coordinates of detected peaks (row, column) in the spectrogram.
/// * `sample_rate` - The sample rate of the audio signal.
///
/// # Returns
///
/// The estimated BPM based on the time intervals between detected peaks.
fn calculate_bpm(peaks: &[(usize, usize)], sample_rate: u32) -> f32 {
  // Calculate time intervals between peaks
  let intervals: Vec<f32> = peaks.windows(2).map(|w| (w[1].0 - w[0].0) as f32 / sample_rate as f32).collect();

  // Use a histogram to find the most common time interval
  let histogram_bins = 100;
  let histogram = create_histogram(&intervals, histogram_bins);

  // Find the bin with the maximum count
  let max_bin = histogram.iter().enumerate().max_by(|(_, &count1), (_, &count2)| count1.cmp(&count2));

  // Calculate the average interval within the most common bin
  if let Some((bin_index, _)) = max_bin {
    let bin_width = 1.0 / histogram_bins as f32;
    let bin_center = (bin_index as f32 + 0.5) * bin_width;
    let average_interval = bin_center * histogram_bins as f32;

    // Convert average interval to BPM
    if average_interval > 0.0 {
      let bpm = 60.0 / average_interval;
      return bpm;
    }
  }

  // Default return value if calculation fails
  0.0
}
