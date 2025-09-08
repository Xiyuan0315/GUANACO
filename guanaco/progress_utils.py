"""Progress indicators for GUANACO startup"""

import sys
import time
import threading

class Spinner:
    """Animated spinner for long operations"""
    
    def __init__(self, message="Loading"):
        self.message = message
        self.running = False
        self.thread = None
        self.spinner_chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        self.current = 0
        
    def __enter__(self):
        self.start()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()
        
    def start(self):
        """Start the spinner"""
        self.running = True
        self.thread = threading.Thread(target=self._spin)
        self.thread.daemon = True
        self.thread.start()
        
    def _spin(self):
        """Spin animation"""
        while self.running:
            char = self.spinner_chars[self.current % len(self.spinner_chars)]
            sys.stdout.write(f'\r{char} {self.message}')
            sys.stdout.flush()
            self.current += 1
            time.sleep(0.1)
            
    def stop(self, success=True):
        """Stop the spinner"""
        self.running = False
        if self.thread:
            self.thread.join()
        
        # Clear the line and show result
        sys.stdout.write('\r' + ' ' * (len(self.message) + 4) + '\r')
        if success:
            print(f"✓ {self.message}")
        else:
            print(f"✗ {self.message}")
        sys.stdout.flush()


def progress_bar(current, total, message="", width=40):
    """Simple progress bar for known quantities"""
    percent = current / total
    filled = int(width * percent)
    bar = "█" * filled + "░" * (width - filled)
    sys.stdout.write(f'\r{message} [{bar}] {percent*100:.0f}%')
    sys.stdout.flush()
    if current >= total:
        print()  # New line when complete