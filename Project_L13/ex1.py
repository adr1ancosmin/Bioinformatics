import numpy as np


def predict_steps(transition_matrix, initial_vector, num_steps=5):
    """
    Make predictions over discrete time steps using a transition matrix.
    
    Parameters:
    -----------
    transition_matrix : numpy.ndarray
        Square matrix representing transitions between states
    initial_vector : numpy.ndarray
        Initial state vector
    num_steps : int
        Number of discrete steps to predict (default: 5)
    
    Returns:
    --------
    list : List of state vectors for each step (including initial state)
    """
    # Validate inputs
    if transition_matrix.shape[0] != transition_matrix.shape[1]:
        raise ValueError("Transition matrix must be square")
    
    if len(initial_vector) != transition_matrix.shape[0]:
        raise ValueError("Initial vector size must match matrix dimensions")
    
    # Store predictions
    predictions = []
    current_state = initial_vector.copy()
    
    # Store initial state (Day 0)
    predictions.append(current_state.copy())
    
    # Calculate predictions for each step
    for step in range(1, num_steps + 1):
        current_state = transition_matrix @ current_state
        predictions.append(current_state.copy())
    
    return predictions


def display_predictions(predictions, step_name="Day"):
    """
    Display predictions in a readable format.
    
    Parameters:
    -----------
    predictions : list
        List of state vectors
    step_name : str
        Name for each step (e.g., "Day", "Step", "Time")
    """
    print(f"\n{'='*50}")
    print(f"Predictions over {len(predictions)-1} {step_name}s")
    print(f"{'='*50}\n")
    
    for i, state in enumerate(predictions):
        print(f"{step_name} {i}:")
        print(f"  State vector: {state}")
        print(f"  Total: {np.sum(state):.4f}")
        print()


def main():
    # Example 1: Simple 2x2 transition matrix
    print("Example 1: 2-state system (e.g., Weather prediction)")
    print("-" * 50)
    
    # Transition matrix (rows sum to 1 for probability interpretation)
    # State 0 -> State 0: 0.7, State 0 -> State 1: 0.3
    # State 1 -> State 0: 0.4, State 1 -> State 1: 0.6
    transition_matrix = np.array([
        [0.7, 0.3],
        [0.4, 0.6]
    ])
    
    # Initial state vector (e.g., 100 units in state 0, 0 in state 1)
    initial_vector = np.array([100.0, 0.0])
    
    # Make predictions
    predictions = predict_steps(transition_matrix, initial_vector, num_steps=5)
    display_predictions(predictions, step_name="Day")
    
    
    # Example 2: 3x3 transition matrix
    print("\n" + "="*50)
    print("Example 2: 3-state system (e.g., Population dynamics)")
    print("-" * 50)
    
    # 3x3 transition matrix
    transition_matrix_3x3 = np.array([
        [0.5, 0.3, 0.2],
        [0.2, 0.6, 0.2],
        [0.1, 0.3, 0.6]
    ])
    
    # Initial state vector
    initial_vector_3x3 = np.array([50.0, 30.0, 20.0])
    
    # Make predictions
    predictions_3x3 = predict_steps(transition_matrix_3x3, initial_vector_3x3, num_steps=5)
    display_predictions(predictions_3x3, step_name="Day")


def custom_prediction():
    """
    Function for custom user input.
    Uncomment and modify this section to use your own matrix and vector.
    """
    print("\n" + "="*50)
    print("Custom Prediction")
    print("-" * 50)
    
    # Define your custom transition matrix here
    # Example: 2x2 matrix
    custom_matrix = np.array([
        [0.8, 0.2],
        [0.3, 0.7]
    ])
    
    # Define your custom initial vector here
    custom_initial = np.array([100.0, 50.0])
    
    # Make predictions
    custom_predictions = predict_steps(custom_matrix, custom_initial, num_steps=5)
    display_predictions(custom_predictions, step_name="Day")


if __name__ == "__main__":
    # Run examples
    main()
    
    # Uncomment the line below to run your custom prediction
    # custom_prediction()
