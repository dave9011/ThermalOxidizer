
package main;

import javax.swing.SwingUtilities;

public class Main {
    
    public static void main(String[] args) {
        
        SwingUtilities.invokeLater( new Runnable(){

            @Override
            public void run() {
                
                try{
                    Thread.sleep(3000);
                } catch(Exception ex){
                }
                
                final ApplicationGUI frame = new ApplicationGUI();
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);               
                
            }
         
        });        
        
    }
    
}
