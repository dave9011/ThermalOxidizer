
package main;

import javax.swing.SwingUtilities;

public class Main {
    
    public static void main(String[] args) {
        
        SwingUtilities.invokeLater( new Runnable(){

            @Override
            public void run() {
                
                try{
                    Thread.sleep(3000);  //2.7sec
                } catch(Exception ex){
                }
                
                final ApplicationGUI frame = new ApplicationGUI();
                //frame.setSize(1400, 650);
                //frame.setExtendedState(frame.getExtendedState()|JFrame.MAXIMIZED_BOTH );
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);               
                
            }
         
        });        
        
    }
    
}
