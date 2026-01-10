"""
Dash Application for 3D Signal Discriminator
=============================================

Interactive web application for signal discrimination
using Plotly Dash.
"""

import numpy as np
from typing import Optional, Dict, Any
import json


def create_dash_app(
    spectrogram_data: Optional[Dict[str, Any]] = None,
    host: str = '127.0.0.1',
    port: int = 8050,
    debug: bool = False
):
    """
    Create Dash application for signal discrimination.

    Args:
        spectrogram_data: Initial spectrogram data (optional)
        host: Server host
        port: Server port
        debug: Enable debug mode

    Returns:
        Dash app instance
    """
    try:
        import dash
        from dash import dcc, html, Input, Output, State
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError(
            "Dash is required for the web application. "
            "Install with: pip install dash plotly"
        )

    # Create app
    app = dash.Dash(
        __name__,
        title='VMS Signal Discriminator',
        suppress_callback_exceptions=True
    )

    # Layout
    app.layout = html.Div([
        # Header
        html.Div([
            html.H1('VMS Signal Discriminator', className='header-title'),
            html.P('Interactive 3D signal visualization and filtering'),
        ], className='header'),

        # Main content
        html.Div([
            # Left panel - Controls
            html.Div([
                html.H3('Controls'),

                # File upload
                html.Div([
                    html.Label('Load Audio File'),
                    dcc.Upload(
                        id='upload-audio',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select File')
                        ]),
                        style={
                            'width': '100%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px 0'
                        },
                    ),
                ]),

                # Threshold slider
                html.Div([
                    html.Label('Signal Threshold (dB)'),
                    dcc.Slider(
                        id='threshold-slider',
                        min=-80,
                        max=-10,
                        step=1,
                        value=-40,
                        marks={i: str(i) for i in range(-80, -9, 10)},
                    ),
                ], style={'margin': '20px 0'}),

                # Frequency range
                html.Div([
                    html.Label('Frequency Range (Hz)'),
                    dcc.RangeSlider(
                        id='freq-slider',
                        min=0,
                        max=22050,
                        step=100,
                        value=[80, 8000],
                        marks={i: str(i) for i in [0, 1000, 5000, 10000, 20000]},
                    ),
                ], style={'margin': '20px 0'}),

                # Presets
                html.Div([
                    html.Label('Filter Presets'),
                    dcc.Dropdown(
                        id='preset-dropdown',
                        options=[
                            {'label': 'Voice Isolation', 'value': 'voice'},
                            {'label': 'Music', 'value': 'music'},
                            {'label': 'Remove Hum', 'value': 'hum'},
                            {'label': 'Custom', 'value': 'custom'},
                        ],
                        value='voice'
                    ),
                ], style={'margin': '20px 0'}),

                # Actions
                html.Div([
                    html.Button('Detect Blocks', id='detect-btn', n_clicks=0,
                               style={'margin': '5px'}),
                    html.Button('Apply Filter', id='apply-btn', n_clicks=0,
                               style={'margin': '5px'}),
                    html.Button('Export Audio', id='export-btn', n_clicks=0,
                               style={'margin': '5px'}),
                ], style={'margin': '20px 0'}),

                # Block list
                html.Div([
                    html.Label('Detected Blocks'),
                    html.Div(id='block-list', children=[]),
                ], style={'margin': '20px 0'}),

            ], className='control-panel', style={'width': '300px', 'padding': '20px'}),

            # Right panel - Visualization
            html.Div([
                # 3D View
                dcc.Graph(
                    id='spectrogram-3d',
                    style={'height': '500px'},
                    config={'responsive': True}
                ),

                # 2D Views
                html.Div([
                    dcc.Graph(
                        id='spectrogram-2d',
                        style={'height': '300px', 'width': '50%', 'display': 'inline-block'}
                    ),
                    dcc.Graph(
                        id='spectrum-avg',
                        style={'height': '300px', 'width': '50%', 'display': 'inline-block'}
                    ),
                ]),

            ], className='visualization-panel', style={'flex': '1', 'padding': '20px'}),

        ], style={'display': 'flex'}),

        # Hidden stores
        dcc.Store(id='spectrogram-store'),
        dcc.Store(id='blocks-store'),
        dcc.Store(id='mask-store'),

        # Download component
        dcc.Download(id='download-audio'),

    ], className='app-container')

    # Callbacks
    @app.callback(
        Output('spectrogram-store', 'data'),
        Input('upload-audio', 'contents'),
        State('upload-audio', 'filename'),
        prevent_initial_call=True
    )
    def load_audio(contents, filename):
        """Load and process uploaded audio."""
        if contents is None:
            return None

        import base64
        import io
        from scipy.io import wavfile

        # Decode base64
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)

        # Read audio
        audio_io = io.BytesIO(decoded)
        try:
            sr, audio = wavfile.read(audio_io)
            if audio.dtype == np.int16:
                audio = audio.astype(np.float32) / 32768.0
        except Exception as e:
            return {'error': str(e)}

        # Convert to mono
        if len(audio.shape) > 1:
            audio = np.mean(audio, axis=1)

        # Compute spectrogram
        from .visualizer import create_3d_spectrogram
        spec = create_3d_spectrogram(audio, sample_rate=sr)

        return {
            'magnitude': spec.magnitude.tolist(),
            'frequencies': spec.frequencies.tolist(),
            'times': spec.times.tolist(),
            'sample_rate': sr,
            'audio': audio.tolist()
        }

    @app.callback(
        Output('spectrogram-3d', 'figure'),
        Output('spectrogram-2d', 'figure'),
        Output('spectrum-avg', 'figure'),
        Input('spectrogram-store', 'data'),
        Input('threshold-slider', 'value'),
        Input('freq-slider', 'value'),
        Input('blocks-store', 'data'),
    )
    def update_visualization(spec_data, threshold, freq_range, blocks_data):
        """Update all visualizations."""
        # Default empty figures
        empty_fig = go.Figure()
        empty_fig.update_layout(
            title='Load an audio file to begin',
            height=400
        )

        if spec_data is None or 'error' in spec_data:
            return empty_fig, empty_fig, empty_fig

        magnitude = np.array(spec_data['magnitude'])
        frequencies = np.array(spec_data['frequencies'])
        times = np.array(spec_data['times'])

        # Frequency mask
        freq_mask = (frequencies >= freq_range[0]) & (frequencies <= freq_range[1])
        freq_display = frequencies[freq_mask]
        mag_display = magnitude[freq_mask, :]

        # 3D Surface
        T, F = np.meshgrid(times, freq_display)
        fig_3d = go.Figure(data=[
            go.Surface(
                x=T, y=F, z=mag_display,
                colorscale='Viridis',
                opacity=0.8
            )
        ])
        fig_3d.update_layout(
            title='3D Spectrogram',
            scene=dict(
                xaxis_title='Time (s)',
                yaxis_title='Frequency (Hz)',
                zaxis_title='Amplitude (dB)',
            ),
            height=500
        )

        # Add threshold plane
        fig_3d.add_trace(go.Surface(
            x=T, y=F,
            z=np.full_like(mag_display, threshold),
            opacity=0.3,
            colorscale=[[0, 'red'], [1, 'red']],
            showscale=False,
            name='Threshold'
        ))

        # 2D Heatmap
        fig_2d = go.Figure(data=[
            go.Heatmap(
                x=times,
                y=freq_display,
                z=mag_display,
                colorscale='Viridis'
            )
        ])
        fig_2d.update_layout(
            title='Time-Frequency View',
            xaxis_title='Time (s)',
            yaxis_title='Frequency (Hz)',
            height=300
        )

        # Average spectrum
        avg_spec = mag_display.mean(axis=1)
        fig_spec = go.Figure(data=[
            go.Scatter(x=freq_display, y=avg_spec, mode='lines', name='Average'),
        ])

        # Add threshold line
        fig_spec.add_hline(y=threshold, line_dash='dash', line_color='red',
                          annotation_text='Threshold')

        fig_spec.update_layout(
            title='Average Spectrum',
            xaxis_title='Frequency (Hz)',
            yaxis_title='Amplitude (dB)',
            height=300
        )

        return fig_3d, fig_2d, fig_spec

    @app.callback(
        Output('blocks-store', 'data'),
        Output('block-list', 'children'),
        Input('detect-btn', 'n_clicks'),
        State('spectrogram-store', 'data'),
        State('threshold-slider', 'value'),
        State('freq-slider', 'value'),
        prevent_initial_call=True
    )
    def detect_blocks(n_clicks, spec_data, threshold, freq_range):
        """Detect signal blocks."""
        if spec_data is None:
            return None, []

        from .visualizer import SignalDiscriminator, Spectrogram3D

        # Reconstruct spectrogram object
        spec = Spectrogram3D(
            magnitude=np.array(spec_data['magnitude']),
            frequencies=np.array(spec_data['frequencies']),
            times=np.array(spec_data['times']),
            sample_rate=spec_data['sample_rate'],
            frame_size=2048,
            hop_size=512
        )

        discriminator = SignalDiscriminator(spec)
        blocks = discriminator.detect_signal_blocks(
            threshold_db=threshold,
            min_freq=freq_range[0],
            max_freq=freq_range[1]
        )

        # Create block list items
        block_items = []
        for block in blocks[:20]:  # Limit display
            block_items.append(
                html.Div([
                    html.Strong(f"Block {block['id']}: "),
                    html.Span(
                        f"{block['freq_range'][0]:.0f}-{block['freq_range'][1]:.0f} Hz, "
                        f"{block['duration']:.2f}s, "
                        f"{block['peak_amplitude']:.1f} dB"
                    )
                ], style={'padding': '5px', 'borderBottom': '1px solid #ccc'})
            )

        # Convert blocks for storage (without numpy arrays)
        blocks_serializable = []
        for b in blocks:
            blocks_serializable.append({
                'id': b['id'],
                'freq_range': list(b['freq_range']),
                'time_range': list(b['time_range']),
                'peak_amplitude': float(b['peak_amplitude']),
                'duration': float(b['duration'])
            })

        return blocks_serializable, block_items

    return app


def run_discriminator(
    audio_path: str = None,
    host: str = '127.0.0.1',
    port: int = 8050,
    debug: bool = False
):
    """
    Run the discriminator web application.

    Args:
        audio_path: Optional path to audio file to load initially
        host: Server host
        port: Server port
        debug: Enable debug mode
    """
    app = create_dash_app(host=host, port=port, debug=debug)

    print(f"\n{'='*50}")
    print("VMS Signal Discriminator")
    print(f"{'='*50}")
    print(f"Open browser at: http://{host}:{port}")
    print("Press Ctrl+C to stop")
    print(f"{'='*50}\n")

    app.run_server(host=host, port=port, debug=debug)


if __name__ == '__main__':
    run_discriminator(debug=True)
